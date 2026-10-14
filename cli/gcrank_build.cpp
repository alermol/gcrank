#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

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
    "                      never by extension)\n"
    "  --out <dir>         Directory that receives <NAME>.gcrank/\n"
    "\n"
    "Options:\n"
    "  --name <name>       Index name. Default: the input file name with a\n"
    "                      compression suffix and a .fa/.fasta/.fna suffix removed\n"
    "  --chunk-size <N>    Bases per parallel work unit. Default: 100000000\n"
    "  --threads <N>       Worker threads. Default: min(8, CPU count)\n"
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
  uint64_t chunk = 100000000ull;
  unsigned threads = gcrank::default_threads();
  bool overwrite = false;
  bool verbose = false;
  bool help = false;
  bool version = false;
};

uint64_t parse_uint(const std::string& option, const std::string& text, uint64_t maximum) {
  if (text.empty()) throw gcrank::UsageError(option + " requires a value");
  uint64_t value = 0;
  for (const char c : text) {
    if (c < '0' || c > '9')
      throw gcrank::UsageError(option + ": not a non-negative number: '" + text + "'");
    const uint64_t digit = static_cast<uint64_t>(c - '0');
    if (value > (maximum - digit) / 10)
      throw gcrank::UsageError(option + ": value out of range: '" + text + "'");
    value = value * 10 + digit;
  }
  return value;
}

Options parse_args(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    std::string value;
    bool has_value = false;
    if (arg.rfind("--", 0) == 0) {
      const std::size_t eq = arg.find('=');
      if (eq != std::string::npos) {
        value = arg.substr(eq + 1);
        arg = arg.substr(0, eq);
        has_value = true;
      }
    }
    auto take_value = [&](const std::string& option) -> std::string {
      if (has_value) return value;
      if (i + 1 >= argc) throw gcrank::UsageError(option + " requires a value");
      return std::string(argv[++i]);
    };
    auto no_value = [&](const std::string& option) {
      if (has_value) throw gcrank::UsageError(option + " does not take a value");
    };

    if (arg == "--fasta") {
      options.fasta = take_value(arg);
    } else if (arg == "--out") {
      options.out = take_value(arg);
    } else if (arg == "--name") {
      options.name = take_value(arg);
    } else if (arg == "--chunk-size") {
      options.chunk = parse_uint(arg, take_value(arg), 1000000000000ull);
    } else if (arg == "--threads") {
      const uint64_t threads = parse_uint(arg, take_value(arg), 1000000ull);
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
      throw gcrank::UsageError("unknown option '" + arg + "'");
    }
  }
  if (!options.help && !options.version) {
    if (options.fasta.empty()) throw gcrank::UsageError("--fasta is required");
    if (options.out.empty()) throw gcrank::UsageError("--out is required");
  }
  return options;
}

bool ends_with_ci(const std::string& text, const std::string& suffix) {
  if (text.size() < suffix.size()) return false;
  const std::size_t offset = text.size() - suffix.size();
  for (std::size_t i = 0; i < suffix.size(); ++i) {
    char a = text[offset + i];
    char b = suffix[i];
    if (a >= 'A' && a <= 'Z') a = static_cast<char>(a + 32);
    if (b >= 'A' && b <= 'Z') b = static_cast<char>(b + 32);
    if (a != b) return false;
  }
  return true;
}

// "reads.fa.gz" -> "reads", "genome" -> "genome".
std::string default_index_name(const std::string& fasta_path) {
  std::string base = std::filesystem::path(fasta_path).filename().string();
  for (const char* suffix : {".gz", ".gzip", ".bgz", ".bgzf"}) {
    if (ends_with_ci(base, suffix)) {
      base.resize(base.size() - std::strlen(suffix));
      break;
    }
  }
  for (const char* suffix : {".fa", ".fasta", ".fna", ".fas", ".fsa", ".txt"}) {
    if (ends_with_ci(base, suffix)) {
      base.resize(base.size() - std::strlen(suffix));
      break;
    }
  }
  if (base.empty()) base = "index";
  return base;
}

}  // namespace

int main(int argc, char** argv) {
  Options options;
  try {
    options = parse_args(argc, argv);
  } catch (const gcrank::UsageError& error) {
    std::fprintf(stderr, "gcrank-build: error: %s\n", error.what());
    std::fprintf(stderr, "Try 'gcrank-build --help' for more information.\n");
    return 1;
  }

  if (options.help) {
    std::fputs(kUsage, stdout);
    return 0;
  }
  if (options.version) {
    std::printf("gcrank-build %s (format %s)\n", gcrank::kBuilderVersion, gcrank::kFormatVersion);
    return 0;
  }

  try {
    const std::string name =
        options.name.empty() ? default_index_name(options.fasta) : options.name;
    gcrank::BuildOptions build_options;
    build_options.chunk_bits = options.chunk;
    build_options.threads = options.threads;
    build_options.verbose = options.verbose;

    gcrank::GcRank index;
    index.build(options.fasta, build_options);
    index.save(options.out, name, options.overwrite);

    if (options.verbose) {
      const std::filesystem::path dir = std::filesystem::path(options.out) / (name + ".gcrank");
      uint64_t bytes = 0;
      std::error_code ec;
      for (const auto& entry : std::filesystem::recursive_directory_iterator(dir, ec))
        if (entry.is_regular_file()) bytes += entry.file_size(ec);
      std::fprintf(stderr, "[gcrank] wrote %s (%llu bytes)\n", dir.string().c_str(),
                   static_cast<unsigned long long>(bytes));
      std::fprintf(stderr, "[gcrank] %llu bases, %llu contigs, GC fraction %.6f\n",
                   static_cast<unsigned long long>(index.total_bases()),
                   static_cast<unsigned long long>(index.contigs().size()),
                   index.gc_fraction(0, index.total_bases(), false));
    }
    return 0;
  } catch (const gcrank::UsageError& error) {
    std::fprintf(stderr, "gcrank-build: error: %s\n", error.what());
    return 1;
  } catch (const gcrank::InputError& error) {
    std::fprintf(stderr, "gcrank-build: input error: %s\n", error.what());
    return 2;
  } catch (const gcrank::OutputError& error) {
    std::fprintf(stderr, "gcrank-build: output error: %s\n", error.what());
    return 3;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "gcrank-build: error: %s\n", error.what());
    return 4;
  }
}
