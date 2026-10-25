#include <algorithm>
#include <charconv>
#include <cstdio>
#include <format>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <string_view>

#include "gcrank/detail/print.h"
#include "gcrank/error.h"
#include "gcrank/gc_rank.h"

namespace {

const char* const kUsage =
    "gcrank-query - answer GC-content queries from a gcrank index\n"
    "\n"
    "Usage:\n"
    "  gcrank-query --index <dir> --bed <file> [options]\n"
    "  gcrank-query --index <dir> --window <N> [--step <S>] [options]\n"
    "\n"
    "Required:\n"
    "  --index <dir>       The <NAME>.gcrank directory to query\n"
    "\n"
    "Windows:\n"
    "  --bed <file>        BED intervals (chrom, start, end; 0-based half-open,\n"
    "                      like bedtools and gcrank). \"-\" reads stdin. Blank\n"
    "                      lines and # / track / browser lines are skipped.\n"
    "  --window <N>        Tile every contig with windows of N bases\n"
    "  --step <S>          Tiling step (default: the window size)\n"
    "\n"
    "Options:\n"
    "  --remove-ambiguous  pct_gc denominator is A+C+G+T instead of the window\n"
    "                      length (bedtools nuc has no such column); a window\n"
    "                      with no A/C/G/T base prints nan\n"
    "  --mmap              Map the index files instead of reading them\n"
    "  --no-validate-crc   Skip the checksum pass over the index (faster start;\n"
    "                      sizes, sidecars and rank totals are still checked)\n"
    "  --version           Print version information and exit\n"
    "  --help, -h          Print this help and exit\n"
    "\n"
    "Output (tab separated, one row per window after a # header):\n"
    "  chrom  start  end  pct_gc  num_gc  num_valid  seq_len\n"
    "pct_gc is printed with %.6f and is bit-identical to bedtools nuc's pct_gc\n"
    "with the default denominator.\n"
    "\n"
    "Exit codes: 0 success, 1 usage, 2 input, 3 output, 4 internal\n";

struct Options {
  std::string index;
  std::string bed;
  uint64_t window = 0;
  uint64_t step = 0;
  bool remove_ambiguous = false;
  bool mmap = false;
  bool validate_crc = true;
  bool help = false;
  bool version = false;
};

uint64_t parse_uint(std::string_view option, std::string_view text) {
  if (text.empty()) throw gcrank::UsageError(std::string(option) + " requires a value");
  uint64_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw gcrank::UsageError(std::string(option) + ": not a non-negative number: '" +
                             std::string(text) + "'");
  return value;
}

Options parse_args(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    std::string_view arg = argv[i];
    std::string_view value;
    bool has_value = false;
    if (arg.starts_with("--")) {
      const std::size_t eq = arg.find('=');
      if (eq != std::string_view::npos) {
        value = arg.substr(eq + 1);
        arg = arg.substr(0, eq);
        has_value = true;
      }
    }
    auto take_value = [&](std::string_view option) -> std::string {
      if (has_value) return std::string(value);
      if (i + 1 >= argc) throw gcrank::UsageError(std::string(option) + " requires a value");
      return std::string(argv[++i]);
    };
    auto no_value = [&](std::string_view option) {
      if (has_value) throw gcrank::UsageError(std::string(option) + " does not take a value");
    };

    if (arg == "--index") {
      options.index = take_value(arg);
    } else if (arg == "--bed") {
      options.bed = take_value(arg);
    } else if (arg == "--window") {
      options.window = parse_uint(arg, take_value(arg));
    } else if (arg == "--step") {
      options.step = parse_uint(arg, take_value(arg));
    } else if (arg == "--remove-ambiguous") {
      no_value(arg);
      options.remove_ambiguous = true;
    } else if (arg == "--mmap") {
      no_value(arg);
      options.mmap = true;
    } else if (arg == "--no-validate-crc") {
      no_value(arg);
      options.validate_crc = false;
    } else if (arg == "--version") {
      no_value(arg);
      options.version = true;
    } else if (arg == "--help" || arg == "-h") {
      no_value(arg);
      options.help = true;
    } else {
      throw gcrank::UsageError("unknown option '" + std::string(arg) + "'");
    }
  }
  if (!options.help && !options.version) {
    if (options.index.empty()) throw gcrank::UsageError("--index is required");
    if (options.bed.empty() && options.window == 0)
      throw gcrank::UsageError("pass either --bed or --window");
    if (!options.bed.empty() && options.window != 0)
      throw gcrank::UsageError("--bed and --window are mutually exclusive");
    if (options.window == 0 && options.step != 0)
      throw gcrank::UsageError("--step requires --window");
    if (options.step != 0 && options.window != 0 && options.step > options.window)
      throw gcrank::UsageError("--step must not exceed --window");
  }
  return options;
}

void write_row(std::FILE* out, std::string_view chrom, uint64_t start, uint64_t end,
               const char* percent, uint64_t gc, uint64_t valid) {
  const std::string line =
      std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", chrom, start, end, percent, gc, valid,
                  end - start);
  if (std::fwrite(line.data(), 1, line.size(), out) != line.size())
    throw gcrank::OutputError("write error on stdout");
}

// GC fraction as a string, or "nan" where the library would reject the window
// (zero length, or no A/C/G/T base with remove_ambiguous).
std::string percent_of(const gcrank::GcRank& index, const gcrank::ContigInfo& contig,
                       uint64_t start, uint64_t end, uint64_t valid, bool remove_ambiguous) {
  if (end == start) return "nan";
  if (remove_ambiguous && valid == 0) return "nan";
  return std::format("{:.6f}",
                     index.gc_fraction(contig.name, start, end, remove_ambiguous));
}

void query_interval(const gcrank::GcRank& index, std::FILE* out, const std::string& chrom,
                    uint64_t start, uint64_t end, uint64_t line_number, bool remove_ambiguous) {
  const gcrank::ContigInfo* contig = index.contig(chrom);
  if (contig == nullptr)
    throw gcrank::InputError("line " + std::to_string(line_number) + ": unknown contig '" +
                             chrom + "'");
  if (start > end || end > contig->length)
    throw gcrank::InputError("line " + std::to_string(line_number) + ": window [" +
                             std::to_string(start) + ", " + std::to_string(end) +
                             ") leaves contig '" + chrom + "' (length " +
                             std::to_string(contig->length) + ")");
  const uint64_t gc = index.gc_count(chrom, start, end);
  const uint64_t valid = index.valid_count(chrom, start, end);
  const std::string percent = percent_of(index, *contig, start, end, valid, remove_ambiguous);
  write_row(out, chrom, start, end, percent.c_str(), gc, valid);
}

bool is_skippable(std::string_view line) {
  while (!line.empty() && (line.front() == ' ' || line.front() == '\t')) line.remove_prefix(1);
  return line.empty() || line.front() == '#' || line.starts_with("track") ||
         line.starts_with("browser");
}

void query_bed(const gcrank::GcRank& index, std::istream& in, std::FILE* out,
               bool remove_ambiguous) {
  std::string line;
  uint64_t line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (is_skippable(line)) continue;
    std::string_view view(line);
    // Split on any run of tabs or spaces; BED is whitespace separated in
    // practice, and extra columns (name, score, ...) are ignored.
    const auto next_field = [&view]() -> std::string_view {
      std::size_t begin = 0;
      while (begin < view.size() && (view[begin] == '\t' || view[begin] == ' ')) ++begin;
      std::size_t end = begin;
      while (end < view.size() && view[end] != '\t' && view[end] != ' ') ++end;
      const std::string_view field = view.substr(begin, end - begin);
      view.remove_prefix(end);
      return field;
    };
    const std::string_view chrom = next_field();
    const std::string_view start_text = next_field();
    const std::string_view end_text = next_field();
    if (chrom.empty() || start_text.empty() || end_text.empty())
      throw gcrank::InputError("line " + std::to_string(line_number) +
                               ": expected at least '<chrom> <start> <end>'");
    if (start_text == "*") continue;  // bedtools' "any contig" wildcard
    const auto field_number = [&](std::string_view text, const char* what) {
      uint64_t value = 0;
      const std::from_chars_result parsed =
          std::from_chars(text.data(), text.data() + text.size(), value);
      if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw gcrank::InputError("line " + std::to_string(line_number) + ": " + what +
                                 " is not a non-negative number: '" + std::string(text) + "'");
      return value;
    };
    const uint64_t start = field_number(start_text, "start");
    const uint64_t end = field_number(end_text, "end");
    query_interval(index, out, std::string(chrom), start, end, line_number, remove_ambiguous);
  }
  if (in.bad()) throw gcrank::InputError("read error on the BED input");
}

void query_tiling(const gcrank::GcRank& index, uint64_t window, uint64_t step, std::FILE* out,
                  bool remove_ambiguous) {
  for (const gcrank::ContigInfo& contig : index.contigs()) {
    for (uint64_t start = 0; start < contig.length; start += step) {
      const uint64_t end = std::min(start + window, contig.length);
      const uint64_t gc = index.gc_count(contig.name, start, end);
      const uint64_t valid = index.valid_count(contig.name, start, end);
      const std::string percent =
          percent_of(index, contig, start, end, valid, remove_ambiguous);
      write_row(out, contig.name, start, end, percent.c_str(), gc, valid);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  Options options;
  try {
    options = parse_args(argc, argv);
  } catch (const gcrank::UsageError& error) {
    gcrank::detail::print(stderr, "gcrank-query: error: {}\n", error.what());
    gcrank::detail::print(stderr, "Try 'gcrank-query --help' for more information.\n");
    return 1;
  }

  if (options.help) {
    gcrank::detail::print(stdout, "{}", kUsage);
    return 0;
  }
  if (options.version) {
    gcrank::detail::print("gcrank-query {} (format {})\n", gcrank::kBuilderVersion,
                          gcrank::kFormatVersion);
    return 0;
  }

  try {
    gcrank::OpenOptions open;
    open.mmap = options.mmap;
    open.validate_crc = options.validate_crc;
    std::mutex log_mutex;
    open.log = [&log_mutex](gcrank::LogLevel level, std::string_view message) {
      std::lock_guard<std::mutex> lock(log_mutex);
      gcrank::detail::print(stderr, "gcrank-query: {}{}\n",
                            level == gcrank::LogLevel::Warning ? "warning: " : "", message);
    };

    gcrank::GcRank index;
    index.load(options.index, open);

    std::FILE* out = stdout;
    if (std::fputs("#chrom\tstart\tend\tpct_gc\tnum_gc\tnum_valid\tseq_len\n", out) < 0)
      throw gcrank::OutputError("write error on stdout");

    if (!options.bed.empty()) {
      std::ifstream file;
      std::istream* in = &std::cin;
      if (options.bed != "-") {
        file.open(options.bed, std::ios::binary);
        if (!file) throw gcrank::InputError("cannot open '" + options.bed + "' for reading");
        in = &file;
      }
      query_bed(index, *in, out, options.remove_ambiguous);
    } else {
      const uint64_t step = options.step == 0 ? options.window : options.step;
      query_tiling(index, options.window, step, out, options.remove_ambiguous);
    }
    if (std::fflush(out) != 0) throw gcrank::OutputError("write error on stdout");
    return 0;
  } catch (const gcrank::UsageError& error) {
    gcrank::detail::print(stderr, "gcrank-query: error: {}\n", error.what());
    return 1;
  } catch (const gcrank::InputError& error) {
    gcrank::detail::print(stderr, "gcrank-query: input error: {}\n", error.what());
    return 2;
  } catch (const gcrank::FormatError& error) {
    gcrank::detail::print(stderr, "gcrank-query: index error: {}\n", error.what());
    return 2;
  } catch (const gcrank::OutputError& error) {
    gcrank::detail::print(stderr, "gcrank-query: output error: {}\n", error.what());
    return 3;
  } catch (const gcrank::IoError& error) {
    // e.g. a missing or unreadable index file: an input problem for a query.
    gcrank::detail::print(stderr, "gcrank-query: input error: {}\n", error.what());
    return 2;
  } catch (const std::exception& error) {
    gcrank::detail::print(stderr, "gcrank-query: error: {}\n", error.what());
    return 4;
  }
}
