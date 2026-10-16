#include <algorithm>
#include <charconv>
#include <cstdio>
#include <filesystem>
#include <mutex>
#include "gcrank/detail/print.h"
#include <string>
#include <string_view>

#if !defined(_WIN32)
#include <unistd.h>
#endif

#include "gcrank/error.h"
#include "gcrank/format.h"
#include "gcrank/gc_rank.h"

namespace {

const char* const kUsage =
    "gcrank-build - build a gcrank index (O(1) GC-content queries) from a FASTA file\n"
    "\n"
    "Usage:\n"
    "  gcrank-build --fasta <path> --out <dir> [options]\n"
    "\n"
    "Required:\n"
    "  --fasta <path>      Input FASTA file, plain or gzip (sniffed by magic bytes,\n"
    "                      never by extension). \"-\" reads stdin, which also\n"
    "                      accepts pipes and process substitution\n"
    "  --out <dir>         Directory that receives <NAME>.gcrank/\n"
    "\n"
    "Options:\n"
    "  --name <name>       Index name. Default: the input file name with a\n"
    "                      compression or FASTA suffix removed (stdin for -)\n"
    "  --chunk-size <N>    Bases per parallel work unit. Default: 16000000\n"
    "                      (maximum 1073741824)\n"
    "  --threads <N>       Worker threads. Default: min(8, CPU count) (max 1024)\n"
    "  --overwrite         Replace an existing index\n"
    "  --verbose           Progress and timing on stderr\n"
    "  --version           Print version information and exit\n"
    "  --help, -h          Print this help and exit\n"
    "\n"
    "Exit codes: 0 success, 1 usage, 2 input, 3 output, 4 internal\n";

struct Options {
  std::string fasta;
  std::string out;
  std::string name;
  uint64_t chunk = 16000000ull;
  unsigned threads = gcrank::default_threads();
  bool overwrite = false;
  bool verbose = false;
  bool help = false;
  bool version = false;
};

uint64_t parse_uint(std::string_view option, std::string_view text, uint64_t maximum) {
  if (text.empty()) throw gcrank::UsageError(std::string(option) + " requires a value");
  uint64_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec == std::errc::result_out_of_range ||
      (parsed.ec == std::errc{} && value > maximum))
    throw gcrank::UsageError(std::string(option) + ": value out of range: '" +
                             std::string(text) + "'");
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
      if (has_value)
        throw gcrank::UsageError(std::string(option) + " does not take a value");
    };

    if (arg == "--fasta") {
      options.fasta = take_value(arg);
    } else if (arg == "--out") {
      options.out = take_value(arg);
    } else if (arg == "--name") {
      options.name = take_value(arg);
    } else if (arg == "--chunk-size") {
      // 1 GiB per work unit is already 1.25 GiB of buffers per slot; refusing
      // absurd values beats silently trying to allocate a terabyte.
      options.chunk = parse_uint(arg, take_value(arg), 1ull << 30);
    } else if (arg == "--threads") {
      const uint64_t threads = parse_uint(arg, take_value(arg), 1024ull);
      if (threads == 0) throw gcrank::UsageError("--threads must be >= 1");
      options.threads = static_cast<unsigned>(threads);
    } else if (arg == "--overwrite") {
      no_value(arg);
      options.overwrite = true;
    } else if (arg == "--verbose") {
      no_value(arg);
      options.verbose = true;
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
    if (options.fasta.empty()) throw gcrank::UsageError("--fasta is required");
    if (options.out.empty()) throw gcrank::UsageError("--out is required");
  }
  return options;
}

char to_ascii_lower(char c) {
  return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
}

bool ends_with_ci(std::string_view text, std::string_view suffix) {
  return text.size() >= suffix.size() &&
         std::ranges::equal(text.substr(text.size() - suffix.size()), suffix, std::equal_to{},
                            to_ascii_lower, to_ascii_lower);
}

// "reads.fa.gz" -> "reads", "genome" -> "genome", "-" -> "stdin".
std::string default_index_name(const std::string& fasta_path) {
  if (fasta_path == "-") return "stdin";
  std::string base = std::filesystem::path(fasta_path).filename().string();
  for (const std::string_view suffix : {".gz", ".gzip", ".bgz", ".bgzf"}) {
    if (ends_with_ci(base, suffix)) {
      base.resize(base.size() - suffix.size());
      break;
    }
  }
  for (const std::string_view suffix : {".fa", ".fasta", ".fna", ".fas", ".fsa", ".txt"}) {
    if (ends_with_ci(base, suffix)) {
      base.resize(base.size() - suffix.size());
      break;
    }
  }
  if (base.empty()) base = "index";
  return base;
}

// Available physical memory in bytes, or 0 when it cannot be determined.
uint64_t available_ram_bytes() {
#if defined(_WIN32)
  return 0;
#else
  const long pages = ::sysconf(_SC_AVPHYS_PAGES);
  const long page_size = ::sysconf(_SC_PAGESIZE);
  if (pages <= 0 || page_size <= 0) return 0;
  return static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_size);
#endif
}

}  // namespace

int main(int argc, char** argv) {
  Options options;
  try {
    options = parse_args(argc, argv);
  } catch (const gcrank::UsageError& error) {
    gcrank::detail::print(stderr, "gcrank-build: error: {}\n", error.what());
    gcrank::detail::print(stderr, "Try 'gcrank-build --help' for more information.\n");
    return 1;
  }

  if (options.help) {
    gcrank::detail::print(stdout, "{}", kUsage);
    return 0;
  }
  if (options.version) {
    gcrank::detail::print("gcrank-build {} (format {})\n", gcrank::kBuilderVersion, gcrank::kFormatVersion);
    return 0;
  }

  try {
    const std::string name =
        options.name.empty() ? default_index_name(options.fasta) : options.name;
    gcrank::BuildOptions build_options;
    build_options.chunk_bits = options.chunk;
    build_options.threads = options.threads;
    build_options.verbose = options.verbose;
    // The library never writes to stderr; the CLI installs the stderr sink.
    std::mutex log_mutex;
    build_options.log = [&log_mutex](gcrank::LogLevel level, std::string_view message) {
      std::lock_guard<std::mutex> lock(log_mutex);
      gcrank::detail::print(stderr, "[gcrank] {}{}\n",
                 level == gcrank::LogLevel::Warning ? "warning: " : "", message);
      std::fflush(stderr);
    };

    gcrank::GcRank index;
    // Warn about an obviously oversized footprint before allocating it.
    const uint64_t estimated = (options.chunk * (options.threads + 2) * 5) / 4;
    const uint64_t available = available_ram_bytes();
    if (available != 0 && estimated > available)
      gcrank::detail::print(stderr,
                 "gcrank-build: warning: estimated peak memory {} MiB exceeds the {} MiB "
                 "available; lower --chunk-size or --threads\n",
                 estimated >> 20, available >> 20);
    index.build(options.fasta, build_options);
    index.save(options.out, name, options.overwrite);

    if (options.verbose) {
      const std::filesystem::path dir = std::filesystem::path(options.out) / (name + ".gcrank");
      uint64_t bytes = 0;
      std::error_code ec;
      for (const auto& entry : std::filesystem::recursive_directory_iterator(dir, ec))
        if (entry.is_regular_file()) bytes += entry.file_size(ec);
      gcrank::detail::print(stderr, "[gcrank] wrote {} ({} bytes)\n", dir.string(), bytes);
      gcrank::detail::print(stderr, "[gcrank] {} bases, {} contigs, ", index.total_bases(),
                 index.contigs().size());
      if (index.total_bases() == 0) {
        // gc_fraction() rejects an empty window, and an empty FASTA is valid.
        gcrank::detail::print(stderr, "GC fraction n/a (empty genome)\n");
      } else {
        gcrank::detail::print(stderr, "GC fraction {:.6f}\n",
                   index.gc_fraction(0, index.total_bases(), false));
      }
    }
    return 0;
  } catch (const gcrank::UsageError& error) {
    gcrank::detail::print(stderr, "gcrank-build: error: {}\n", error.what());
    return 1;
  } catch (const gcrank::InputError& error) {
    gcrank::detail::print(stderr, "gcrank-build: input error: {}\n", error.what());
    return 2;
  } catch (const gcrank::OutputError& error) {
    gcrank::detail::print(stderr, "gcrank-build: output error: {}\n", error.what());
    return 3;
  } catch (const std::exception& error) {
    gcrank::detail::print(stderr, "gcrank-build: error: {}\n", error.what());
    return 4;
  }
}
