#include <doctest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "gcrank/detail/bitvector.h"
#include "gcrank/detail/crc32.h"
#include "gcrank/detail/input_stream.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

struct Rng {
  uint64_t state;
  uint64_t next() {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return state >> 33;
  }
  bool bit() { return (next() & 1u) != 0; }
};

std::vector<bool> random_bits(uint64_t count, uint64_t seed) {
  Rng rng{seed};
  std::vector<bool> bits(static_cast<std::size_t>(count));
  for (uint64_t i = 0; i < count; ++i) bits[static_cast<std::size_t>(i)] = rng.bit();
  return bits;
}

std::vector<uint64_t> naive_ranks(const std::vector<bool>& bits) {
  std::vector<uint64_t> ranks(bits.size() + 1, 0);
  for (std::size_t i = 0; i < bits.size(); ++i)
    ranks[i + 1] = ranks[i] + (bits[i] ? uint64_t(1) : uint64_t(0));
  return ranks;
}

BitVector make_bitvector(const std::vector<bool>& bits) {
  BitVector bitvector;
  for (const bool bit : bits) bitvector.push_bit(bit);
  bitvector.finalize();
  return bitvector;
}

// Packs a slice of bits into a standalone MSB-first word block.
std::vector<uint64_t> pack_slice(const std::vector<bool>& bits, uint64_t start, uint64_t count) {
  std::vector<uint64_t> words(static_cast<std::size_t>((count + 63) / 64), 0);
  for (uint64_t i = 0; i < count; ++i)
    if (bits[static_cast<std::size_t>(start + i)])
      words[static_cast<std::size_t>(i / 64)] |= uint64_t(1) << (63 - (i % 64));
  return words;
}

// The bits after `size()` in the last word must be zero (they are padding);
// the valid bits live in the *top* of the word.
bool tail_padding_is_zero(const BitVector& bitvector) {
  if (bitvector.num_words() == 0) return true;
  const unsigned rem = static_cast<unsigned>(bitvector.size() % 64);
  if (rem == 0) return true;
  const uint64_t tail = bitvector.word(bitvector.num_words() - 1);
  return (tail & ~mask_high_bits(rem)) == 0;
}

const uint64_t kBoundaryLengths[] = {0, 1, 63, 64, 65, 511, 512, 513, 4095, 4096, 4097};

}  // namespace

TEST_CASE("bitvector: helpers") {
  CHECK(std::popcount(0u) == 0);
  CHECK(std::popcount(1u) == 1);
  CHECK(std::popcount(~uint64_t(0)) == 64);
  CHECK(std::popcount(0x8000000000000001ull) == 2);
  CHECK(mask_high_bits(0) == 0);
  CHECK(mask_high_bits(1) == 0x8000000000000000ull);
  CHECK(mask_high_bits(2) == 0xC000000000000000ull);
  CHECK(mask_high_bits(63) == 0xFFFFFFFFFFFFFFFEull);
  CHECK(mask_high_bits(64) == ~uint64_t(0));
}

TEST_CASE("bitvector: push_bit, size and bit") {
  BitVector bitvector;
  CHECK(bitvector.size() == 0);
  CHECK(bitvector.num_words() == 0);
  CHECK_THROWS_AS(bitvector.bit(0), QueryError);

  const std::vector<bool> bits = random_bits(1000, 42);
  for (std::size_t i = 0; i < bits.size(); ++i) {
    bitvector.push_bit(bits[i]);
    CHECK(bitvector.size() == i + 1);
  }
  for (std::size_t i = 0; i < bits.size(); ++i) CHECK(bitvector.bit(i) == bits[i]);
  CHECK_THROWS_AS(bitvector.bit(bits.size()), QueryError);
  CHECK_THROWS_AS(bitvector.rank1(0), Error);  // not finalized yet
}

TEST_CASE("bitvector: empty vector") {
  BitVector bitvector;
  bitvector.finalize();
  CHECK(bitvector.size() == 0);
  CHECK(bitvector.num_words() == 0);
  CHECK(bitvector.has_rank());
  CHECK(bitvector.rank1(0) == 0);
  CHECK_THROWS_AS(bitvector.rank1(1), QueryError);
}

TEST_CASE("bitvector: rank1 matches a naive count at every position") {
  for (const uint64_t length : kBoundaryLengths) {
    CAPTURE(length);
    const std::vector<bool> bits = random_bits(length, length * 7919 + 11);
    const BitVector bitvector = make_bitvector(bits);
    const std::vector<uint64_t> ranks = naive_ranks(bits);
    CHECK(bitvector.size() == length);
    for (uint64_t position = 0; position <= length; ++position)
      CHECK(bitvector.rank1(position) == ranks[position]);
    CHECK_THROWS_AS(bitvector.rank1(length + 1), QueryError);
  }
}

TEST_CASE("bitvector: rank1 is monotonic") {
  const std::vector<bool> bits = random_bits(4097, 99);
  const BitVector bitvector = make_bitvector(bits);
  uint64_t previous = 0;
  for (uint64_t position = 0; position <= bitvector.size(); position += 7) {
    const uint64_t rank = bitvector.rank1(position);
    CHECK(rank >= previous);
    previous = rank;
  }
  CHECK(bitvector.rank1(bitvector.size()) == naive_ranks(bits)[bits.size()]);
}

TEST_CASE("bitvector: append_bits matches push_bit across misaligned blocks") {
  for (const uint64_t length : kBoundaryLengths) {
    CAPTURE(length);
    const std::vector<bool> bits = random_bits(length, length * 31 + 7);
    const std::vector<uint64_t> ranks = naive_ranks(bits);

    BitVector bitvector;
    Rng rng{length * 17 + 3};
    uint64_t offset = 0;
    while (offset < length) {
      const uint64_t take = std::min<uint64_t>(1 + rng.next() % 200, length - offset);
      const std::vector<uint64_t> block = pack_slice(bits, offset, take);
      bitvector.append_bits(block, take);
      offset += take;
    }
    CHECK(offset == length);
    CHECK(bitvector.size() == length);
    bitvector.finalize();
    for (uint64_t position = 0; position <= length; ++position)
      CHECK(bitvector.rank1(position) == ranks[position]);
    CHECK(tail_padding_is_zero(bitvector));
  }
}

TEST_CASE("bitvector: append_bits of zero bits is a no-op") {
  BitVector bitvector;
  bitvector.push_bit(true);
  bitvector.append_bits({}, 0);
  CHECK(bitvector.size() == 1);
  CHECK(bitvector.word(0) == 0x8000000000000000ull);
}

TEST_CASE("bitvector: append_bits rejects a source shorter than nbits") {
  BitVector bitvector;
  const std::vector<uint64_t> one_word{0};
  CHECK_THROWS_AS(bitvector.append_bits(one_word, 65), UsageError);
  CHECK(bitvector.size() == 0);
}

TEST_CASE("bitvector: padding bits stay zero") {
  const std::vector<bool> bits = random_bits(300, 5);

  BitVector pushed;
  for (const bool bit : bits) pushed.push_bit(bit);
  CHECK(tail_padding_is_zero(pushed));

  BitVector appended;
  const std::vector<uint64_t> words = pack_slice(bits, 0, bits.size());
  appended.append_bits(words, bits.size());
  CHECK(tail_padding_is_zero(appended));

  BitVector halves;
  halves.append_bits(words, 64);
  const std::vector<uint64_t> rest = pack_slice(bits, 64, bits.size() - 64);
  halves.append_bits(rest, bits.size() - 64);
  CHECK(tail_padding_is_zero(halves));
  CHECK(halves.size() == bits.size());

  halves.finalize();
  const BitVector reference = make_bitvector(bits);
  for (uint64_t position = 0; position <= bits.size(); ++position)
    CHECK(halves.rank1(position) == reference.rank1(position));
}

TEST_CASE("bitvector: save/load/map roundtrip") {
  gctest::TempDir tmp("bitvector");
  const fs::path bv_path = tmp.file("genome.bv");
  const fs::path rb_path = tmp.file("genome.rb");

  for (const uint64_t length : kBoundaryLengths) {
    CAPTURE(length);
    const std::vector<bool> bits = random_bits(length, length * 13 + 3);
    const BitVector bitvector = make_bitvector(bits);

    const FileInfo bv_info = bitvector.save_bv(bv_path.string());
    const FileInfo rb_info = bitvector.save_rank(rb_path.string());
    CHECK(bv_info.size == bv_bytes(length));
    CHECK(rb_info.size == rank_bytes(length));
    CHECK(fs::file_size(bv_path) == bv_bytes(length));
    CHECK(fs::file_size(rb_path) == rank_bytes(length));
    CHECK(bv_info.crc32 == crc32_file(bv_path.string()));
    CHECK(rb_info.crc32 == crc32_file(rb_path.string()));

    BitVector loaded;
    loaded.load_bv(bv_path.string(), length);
    loaded.load_rank(rb_path.string(), length, host_is_little_endian());
    CHECK(loaded.size() == length);
    CHECK(!loaded.is_mapped());
    for (uint64_t position = 0; position <= length; ++position)
      CHECK(loaded.rank1(position) == bitvector.rank1(position));
    for (uint64_t i = 0; i < length; i += 37) CHECK(loaded.bit(i) == bitvector.bit(i));

    BitVector mapped;
    mapped.map_bv(MappedFile::map(bv_path.string()), length);
    mapped.map_rank(MappedFile::map(rb_path.string()), length, host_is_little_endian());
    CHECK(mapped.is_mapped());
    CHECK(mapped.size() == length);
    CHECK(mapped.num_words() == bv_bytes(length) / 8);
    for (uint64_t position = 0; position <= length; ++position)
      CHECK(mapped.rank1(position) == bitvector.rank1(position));
    for (uint64_t i = 0; i < length; i += 37) CHECK(mapped.bit(i) == bitvector.bit(i));

    // Copying a mapped vector materializes its data.
    BitVector copy = mapped;
    CHECK(copy.size() == length);
    CHECK(!copy.is_mapped());
    CHECK(copy.rank1(length) == bitvector.rank1(length));
  }
}

TEST_CASE("bitvector: rank arrays have the documented extra entry") {
  const std::vector<bool> bits = random_bits(4096, 8);
  const BitVector bitvector = make_bitvector(bits);
  // 4096 bits = exactly one superblock and eight basic blocks, so both arrays
  // need their trailing sentinel entry.
  CHECK(rank_superblock_entries(4096) == 2);
  CHECK(rank_basic_block_entries(4096) == 9);
  CHECK(bitvector.rank1(4096) == bitvector.rank1(4095) + (bits[4095] ? 1u : 0u));
}

TEST_CASE("bitvector: rank files are little-endian on every host") {
  gctest::TempDir tmp("bitvector-endian-write");
  const std::vector<bool> bits = random_bits(4097, 7);
  const BitVector bitvector = make_bitvector(bits);
  const fs::path rb_path = tmp.file("genome.rb");
  bitvector.save_rank(rb_path.string());

  const std::vector<uint8_t> bytes = gctest::read_bytes(rb_path);
  const uint64_t nsb = rank_superblock_entries(4097);
  const uint64_t nbb = rank_basic_block_entries(4097);
  REQUIRE(bytes.size() == rank_bytes(4097));
  const auto le64 = [&bytes](std::size_t at) {
    uint64_t value = 0;
    for (int b = 0; b < 8; ++b) value |= uint64_t(bytes[at + static_cast<std::size_t>(b)]) << (8 * b);
    return value;
  };
  const auto le16 = [&bytes](std::size_t at) {
    return uint32_t(bytes[at]) | (uint32_t(bytes[at + 1]) << 8);
  };
  // Decode with an explicit little-endian reader (no host assumptions) and
  // compare against the rank arithmetic at every entry, including the padding
  // entries past the end of the vector.
  for (uint64_t i = 0; i < nsb; ++i) {
    const uint64_t pos = std::min<uint64_t>(i * kSuperblockBits, 4097);
    CHECK(le64(static_cast<std::size_t>(8 * i)) == bitvector.rank1(pos));
  }
  for (uint64_t j = 0; j < nbb; ++j) {
    const uint64_t pos = std::min<uint64_t>(j * kBasicBlockBits, 4097);
    const uint64_t superblock = (pos / kSuperblockBits) * kSuperblockBits;
    CHECK(le16(static_cast<std::size_t>(8 * nsb + 2 * j)) ==
          bitvector.rank1(pos) - bitvector.rank1(superblock));
  }
}

TEST_CASE("bitvector: a big-endian rank file loads and maps identically") {
  gctest::TempDir tmp("bitvector-endian-read");
  const std::vector<bool> bits = random_bits(4097, 42);
  const BitVector bitvector = make_bitvector(bits);
  const fs::path rb_path = tmp.file("genome.rb");
  const fs::path bv_path = tmp.file("genome.bv");
  bitvector.save_bv(bv_path.string());
  bitvector.save_rank(rb_path.string());

  // Byte-swap every scalar, as a writer on a big-endian host would store it.
  std::vector<uint8_t> bytes = gctest::read_bytes(rb_path);
  const uint64_t nsb = rank_superblock_entries(4097);
  const uint64_t nbb = rank_basic_block_entries(4097);
  for (uint64_t i = 0; i < nsb; ++i)
    std::reverse(bytes.begin() + static_cast<std::ptrdiff_t>(8 * i),
                 bytes.begin() + static_cast<std::ptrdiff_t>(8 * i + 8));
  for (uint64_t j = 0; j < nbb; ++j)
    std::swap(bytes[static_cast<std::size_t>(8 * nsb + 2 * j)],
              bytes[static_cast<std::size_t>(8 * nsb + 2 * j + 1)]);
  const fs::path big_endian_path = tmp.file("genome_be.rb");
  gctest::write_text(big_endian_path, std::string(bytes.begin(), bytes.end()));

  BitVector loaded;
  loaded.load_bv(bv_path.string(), 4097);
  loaded.load_rank(big_endian_path.string(), 4097, /*stored_little_endian=*/false);
  BitVector mapped;
  mapped.map_bv(MappedFile::map(bv_path.string()), 4097);
  mapped.map_rank(MappedFile::map(big_endian_path.string()), 4097, false);
  for (uint64_t position = 0; position <= 4097; ++position) {
    CHECK(loaded.rank1(position) == bitvector.rank1(position));
    CHECK(mapped.rank1(position) == bitvector.rank1(position));
  }
}

TEST_CASE("bitvector: rank validation rejects a broken basic-block index") {
  gctest::TempDir tmp("bitvector-validation");
  const std::vector<bool> bits = random_bits(4097, 11);
  const BitVector bitvector = make_bitvector(bits);
  const fs::path rb_path = tmp.file("genome.rb");
  bitvector.save_rank(rb_path.string());
  const uint64_t nsb = rank_superblock_entries(4097);
  const uint64_t nbb = rank_basic_block_entries(4097);
  const std::vector<uint8_t> clean = gctest::read_bytes(rb_path);

  // The first entry of every superblock is a rank relative to that superblock
  // and must be zero: entry 8 starts the second superblock (4096 bits).
  REQUIRE(nbb > 8);
  {
    std::vector<uint8_t> bytes = clean;
    const std::size_t at = static_cast<std::size_t>(8 * nsb + 2 * 8);
    bytes[at] = 1;
    bytes[at + 1] = 0;
    gctest::write_text(rb_path, std::string(bytes.begin(), bytes.end()));
    BitVector broken;
    CHECK_THROWS_AS(broken.load_rank(rb_path.string(), 4097, true), FormatError);
  }
  // A decreasing pair inside one superblock is impossible.
  {
    std::vector<uint8_t> bytes = clean;
    // Entry 2 follows entry 1 in the first superblock; force it below entry 1
    // (or wrap to 0xFFFF, which is also out of range).
    const uint16_t previous = static_cast<uint16_t>(bytes[8 * nsb + 2] |
                                                    (bytes[8 * nsb + 3] << 8));
    const uint16_t forbidden = static_cast<uint16_t>(previous - 1);
    bytes[8 * nsb + 4] = static_cast<uint8_t>(forbidden & 0xFFu);
    bytes[8 * nsb + 5] = static_cast<uint8_t>(forbidden >> 8);
    gctest::write_text(rb_path, std::string(bytes.begin(), bytes.end()));
    BitVector broken;
    CHECK_THROWS_AS(broken.load_rank(rb_path.string(), 4097, true), FormatError);
  }
  // The clean file still loads.
  gctest::write_text(rb_path, std::string(clean.begin(), clean.end()));
  BitVector good;
  CHECK_NOTHROW(good.load_rank(rb_path.string(), 4097, true));
}

TEST_CASE("bitvector: corrupt or mismatched files are rejected") {
  gctest::TempDir tmp("bitvector-errors");
  const std::vector<bool> bits = random_bits(513, 21);
  const BitVector bitvector = make_bitvector(bits);
  const fs::path bv_path = tmp.file("genome.bv");
  const fs::path rb_path = tmp.file("genome.rb");
  bitvector.save_bv(bv_path.string());
  bitvector.save_rank(rb_path.string());

  BitVector other;
  // 512 and 513 bits have the same padded size, so use a length that cannot
  // match the file on disk.
  CHECK_THROWS_AS(other.load_bv(bv_path.string(), 4097), FormatError);
  CHECK_THROWS_AS(other.load_rank(rb_path.string(), 4097, true), FormatError);

  // A truncated bitvector file.
  std::vector<uint8_t> bytes = gctest::read_bytes(bv_path);
  bytes.resize(bytes.size() - 8);
  gctest::write_text(bv_path, std::string(bytes.begin(), bytes.end()));
  BitVector truncated;
  CHECK_THROWS_AS(truncated.load_bv(bv_path.string(), 513), FormatError);

  // A rank index whose basic-block entries are impossible.
  gctest::write_text(rb_path, std::string(rank_bytes(513), '\xff'));
  BitVector broken;
  CHECK_THROWS_AS(broken.load_rank(rb_path.string(), 513, true), FormatError);
}
