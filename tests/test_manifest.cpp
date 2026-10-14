#include <doctest.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "gcrank/crc32.h"
#include "gcrank/format.h"
#include "gcrank/manifest.h"
#include "gcrank/sha256.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

Manifest sample_manifest() {
  Manifest manifest;
  manifest.created = "2023-11-14T22:13:20Z";
  manifest.source.path = "tests/data/test_small.fa";
  manifest.source.sha256 = std::string(64, 'a');
  manifest.source.file_bytes = 11530;
  manifest.source.total_bases = 11095;
  manifest.source.num_contigs = 7;
  manifest.endianness = host_endianness_name();
  manifest.files["contigs.tsv"] = FileInfo{128, 0x11223344};
  manifest.files["genome.bv"] = FileInfo{1408, 0xdeadbeef};
  manifest.files["genome.rb"] = FileInfo{250, 0x0badc0de};
  manifest.files["valid.bv"] = FileInfo{1408, 0xfeedface};
  manifest.files["valid.rb"] = FileInfo{250, 0x00c0ffee};
  manifest.semantics = default_semantics();
  return manifest;
}

std::string rename_key(const std::string& text, const std::string& from, const std::string& to) {
  const std::string needle = "\"" + from + "\":";
  const std::string replacement = "\"" + to + "\":";
  const std::size_t position = text.find(needle);
  if (position == std::string::npos) return text;
  std::string result = text;
  result.replace(position, needle.size(), replacement);
  return result;
}

struct KeyLine {
  std::size_t indent = 0;
  std::string key;
};

std::vector<KeyLine> key_lines(const std::string& text) {
  std::vector<KeyLine> keys;
  std::size_t line_begin = 0;
  while (line_begin <= text.size()) {
    std::size_t line_end = text.find('\n', line_begin);
    if (line_end == std::string::npos) line_end = text.size();
    const std::string line = text.substr(line_begin, line_end - line_begin);
    std::size_t indent = 0;
    while (indent < line.size() && line[indent] == ' ') ++indent;
    if (indent < line.size() && line[indent] == '"') {
      const std::size_t closing = line.find("\":", indent + 1);
      if (closing != std::string::npos && closing + 2 < line.size() &&
          line[closing + 2] == ' ') {
        KeyLine entry;
        entry.indent = indent;
        entry.key = line.substr(indent + 1, closing - indent - 1);
        keys.push_back(entry);
      }
    }
    if (line_end == text.size()) break;
    line_begin = line_end + 1;
  }
  return keys;
}

}  // namespace

TEST_CASE("manifest: JSON round trip preserves every field") {
  const Manifest original = sample_manifest();
  const std::string text = original.to_json();
  const Manifest parsed = Manifest::from_json(text, "test");
  CHECK(parsed.format_version == original.format_version);
  CHECK(parsed.builder_version == original.builder_version);
  CHECK(parsed.created == original.created);
  CHECK(parsed.source.path == original.source.path);
  CHECK(parsed.source.sha256 == original.source.sha256);
  CHECK(parsed.source.file_bytes == original.source.file_bytes);
  CHECK(parsed.source.total_bases == original.source.total_bases);
  CHECK(parsed.source.num_contigs == original.source.num_contigs);
  CHECK(parsed.bit_order == original.bit_order);
  CHECK(parsed.endianness == original.endianness);
  CHECK(parsed.superblock_bits == original.superblock_bits);
  CHECK(parsed.basic_block_bits == original.basic_block_bits);
  CHECK(parsed.remove_ambiguous_supported == original.remove_ambiguous_supported);
  CHECK(parsed.files.size() == original.files.size());
  for (const auto& entry : original.files) {
    CHECK(parsed.files.at(entry.first).size == entry.second.size);
    CHECK(parsed.files.at(entry.first).crc32 == entry.second.crc32);
  }
  CHECK(parsed.semantics.at("gc_bases") == "GC");
  CHECK(parsed.semantics.at("valid_bases") == "ACGT");
}

TEST_CASE("manifest: canonical JSON layout") {
  const std::string text = sample_manifest().to_json();
  CHECK(text.back() == '\n');
  CHECK(text.find('\r') == std::string::npos);
  CHECK(text.find("\\u") == std::string::npos);  // no stray escapes
  CHECK(text.rfind('\n') == text.size() - 1);    // exactly one trailing newline

  // Keys are sorted within every object (same indent = same object).
  const std::vector<KeyLine> keys = key_lines(text);
  CHECK(keys.size() >= 20);
  for (std::size_t i = 1; i < keys.size(); ++i)
    if (keys[i].indent == keys[i - 1].indent) CHECK(keys[i].key > keys[i - 1].key);

  // Indentation is a multiple of two spaces, never tabs.
  for (const KeyLine& entry : keys) CHECK(entry.indent % 2 == 0);
  CHECK(text.find('\t') == std::string::npos);

  // Known canonical line for a scalar.
  const std::string version_line =
      "\n  \"format_version\": \"" + std::string(kFormatVersion) + "\",\n";
  CHECK(text.find(version_line) != std::string::npos);
  CHECK(text.find("\n      \"crc32\": \"deadbeef\",\n") != std::string::npos);
}

TEST_CASE("manifest: a different major version is rejected") {
  const std::string current = std::string("\"format_version\": \"") + kFormatVersion + "\"";
  std::string text = sample_manifest().to_json();
  text = rename_key(text, "format_version", "format_version_");  // keep JSON valid
  const std::string renamed = std::string("\"format_version_\": \"") + kFormatVersion + "\"";
  const std::size_t position = text.find(renamed);
  REQUIRE(position != std::string::npos);
  text.replace(position, renamed.size(), "\"format_version\": \"2.0.0\"");
  CHECK_THROWS_AS(Manifest::from_json(text, "test"), FormatError);

  // Same major, newer minor: accepted.
  std::string minor = sample_manifest().to_json();
  const std::size_t at = minor.find(current);
  REQUIRE(at != std::string::npos);
  minor.replace(at, current.size(), "\"format_version\": \"1.9.3\"");
  CHECK(Manifest::from_json(minor, "test").format_version == "1.9.3");

  // Not a version number at all.
  std::string bad = sample_manifest().to_json();
  const std::size_t bad_at = bad.find(current);
  REQUIRE(bad_at != std::string::npos);
  bad.replace(bad_at, current.size(), "\"format_version\": \"one\"");
  CHECK_THROWS_AS(Manifest::from_json(bad, "test"), FormatError);
}

TEST_CASE("manifest: every required key is enforced") {
  const std::string text = sample_manifest().to_json();
  const char* required[] = {
      "format_version", "builder_version", "created", "source",   "index",
      "files",          "semantics",      "bit_order", "endianness",
      "superblock_bits", "basic_block_bits", "remove_ambiguous_supported", "path",
      "sha256",         "total_bases",    "num_contigs", "size", "crc32",
      "gc_bases",       "valid_bases",
  };
  for (const char* key : required) {
    CAPTURE(key);
    const std::string renamed = rename_key(text, key, std::string("x_") + key);
    REQUIRE(renamed != text);
    CHECK_THROWS_AS(Manifest::from_json(renamed, "test"), FormatError);
  }
}

TEST_CASE("manifest: unknown keys are ignored") {
  std::string text = sample_manifest().to_json();
  const std::size_t position = text.find("\"created\":");
  REQUIRE(position != std::string::npos);
  text.insert(position, "\"future_field\": [1, 2, {\"nested\": true}],\n  ");
  const Manifest parsed = Manifest::from_json(text, "test");
  CHECK(parsed.source.total_bases == 11095);
}

TEST_CASE("manifest: unsupported index parameters are rejected") {
  const auto replace_once = [](std::string text, const std::string& from, const std::string& to) {
    const std::size_t position = text.find(from);
    if (position != std::string::npos) text.replace(position, from.size(), to);
    return text;
  };
  const std::string base = sample_manifest().to_json();

  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"bit_order\": \"msb0\"",
                                                   "\"bit_order\": \"lsb0\""), "test"),
                  FormatError);
  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"endianness\": \"" + std::string(host_endianness_name()) + "\"",
                                                   "\"endianness\": \"middle\""), "test"),
                  FormatError);
  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"superblock_bits\": 4096",
                                                   "\"superblock_bits\": 1024"), "test"),
                  FormatError);
  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"remove_ambiguous_supported\": true",
                                                   "\"remove_ambiguous_supported\": false"), "test"),
                  FormatError);
  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"crc32\": \"deadbeef\"",
                                                   "\"crc32\": \"deadbee\""), "test"),
                  FormatError);
  CHECK_THROWS_AS(Manifest::from_json(replace_once(base, "\"size\": 1408",
                                                   "\"size\": \"big\""), "test"),
                  FormatError);
}

TEST_CASE("manifest: malformed JSON is rejected") {
  const std::string text = sample_manifest().to_json();
  CHECK_THROWS_AS(Manifest::from_json("", "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json("{", "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json(text.substr(0, text.size() / 2), "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json(text + "{}", "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json("[1, 2, 3]", "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json("{\"a\": \"\\uD800\"}", "test"), FormatError);
  CHECK_THROWS_AS(Manifest::from_json("{\"a\": 'single'}", "test"), FormatError);
}

TEST_CASE("manifest: trailing whitespace after the document is accepted") {
  const std::string text = sample_manifest().to_json();
  const Manifest parsed = Manifest::from_json(text + "\n", "test");
  CHECK(parsed.source.num_contigs == 7);
}

TEST_CASE("manifest: saving and loading a file") {
  gctest::TempDir tmp("manifest");
  const fs::path path = tmp.file("manifest.json");
  const Manifest original = sample_manifest();
  original.save(path.string());
  const Manifest loaded = Manifest::load(path.string());
  CHECK(loaded.to_json() == original.to_json());
  CHECK(loaded.created == original.created);
  CHECK_THROWS_AS(Manifest::load(tmp.file("missing.json").string()), IoError);
}

TEST_CASE("manifest: created is derived from SOURCE_DATE_EPOCH") {
  {
    const gctest::EnvGuard guard("SOURCE_DATE_EPOCH", "1700000000");
    CHECK(created_from_env() == "2023-11-14T22:13:20Z");
  }
  {
    const gctest::EnvGuard guard("SOURCE_DATE_EPOCH", "0");
    CHECK(created_from_env() == "1970-01-01T00:00:00Z");
  }
  {
    const gctest::EnvGuard guard("SOURCE_DATE_EPOCH", "not-a-number");
    CHECK_THROWS_AS(created_from_env(), UsageError);
  }
}

TEST_CASE("crc32: known answers and file checksums") {
  // The standard check value for CRC-32 (IEEE 802.3).
  CHECK(CRC32::compute("123456789", 9) == 0xCBF43926u);
  CHECK(CRC32::compute("", 0) == 0x00000000u);
  CHECK(crc32_hex(0xCBF43926u) == "cbf43926");
  CHECK(crc32_hex(0u) == "00000000");

  // Streaming updates agree with a single update.
  CRC32 streaming;
  const std::string text = "the quick brown fox jumps over the lazy dog";
  streaming.update(text.data(), 5);
  streaming.update(text.data() + 5, text.size() - 5);
  CHECK(streaming.value() == CRC32::compute(text.data(), text.size()));

  gctest::TempDir tmp("crc32");
  const fs::path path = tmp.file("data.bin");
  gctest::write_text(path, text);
  CHECK(crc32_file(path.string()) == CRC32::compute(text.data(), text.size()));
  CHECK_THROWS_AS(crc32_file(tmp.file("nope.bin").string()), IoError);
}

TEST_CASE("sha256: known answers, streaming and files") {
  CHECK(sha256_hex("", 0) ==
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  CHECK(sha256_hex("abc", 3) ==
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  const char* fifty_six = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  CHECK(sha256_hex(fifty_six, 56) ==
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
  const std::string million(1000000, 'a');
  CHECK(sha256_hex(million.data(), million.size()) ==
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");

  // Chunked updates cross the 64-byte block boundary in many places.
  SHA256 streaming;
  std::size_t offset = 0;
  std::size_t step = 1;
  while (offset < million.size()) {
    const std::size_t take = std::min(step, million.size() - offset);
    streaming.update(million.data() + offset, take);
    offset += take;
    step = step * 3 + 1;
  }
  uint8_t digest[32];
  streaming.finish(digest);
  std::string hex;
  for (int i = 0; i < 32; ++i) {
    static const char* digits = "0123456789abcdef";
    hex.push_back(digits[digest[i] >> 4]);
    hex.push_back(digits[digest[i] & 0x0F]);
  }
  CHECK(hex == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");

  gctest::TempDir tmp("sha256");
  const fs::path path = tmp.file("data.bin");
  gctest::write_text(path, million);
  CHECK(sha256_file_hex(path.string()) ==
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
  CHECK_THROWS_AS(sha256_file_hex(tmp.file("nope.bin").string()), IoError);
}
