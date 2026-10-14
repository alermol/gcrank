#include "gcrank/bitvector.h"

#include <algorithm>
#include <cstring>
#include <fstream>

#include "gcrank/crc32.h"

namespace gcrank {

namespace {
inline uint64_t load_le64(const uint8_t* p) {
  return static_cast<uint64_t>(p[0]) | (static_cast<uint64_t>(p[1]) << 8) |
         (static_cast<uint64_t>(p[2]) << 16) | (static_cast<uint64_t>(p[3]) << 24) |
         (static_cast<uint64_t>(p[4]) << 32) | (static_cast<uint64_t>(p[5]) << 40) |
         (static_cast<uint64_t>(p[6]) << 48) | (static_cast<uint64_t>(p[7]) << 56);
}
inline void store_le64(uint8_t* p, uint64_t v) {
  for (int i = 0; i < 8; ++i) p[i] = static_cast<uint8_t>(v >> (8 * i));
}
inline uint16_t bswap16(uint16_t v) { return static_cast<uint16_t>((v >> 8) | (v << 8)); }
}  // namespace

// ---- build -----------------------------------------------------------------

void BitVector::clear() {
  words_.clear();
  words_view_ = nullptr;
  num_words_ = 0;
  size_ = 0;
  clear_rank();
  maps_.clear();
  mapped_ = false;
}

void BitVector::clear_rank() {
  sb_owned_.clear();
  bb_owned_.clear();
  sb_view_ = nullptr;
  bb_view_ = nullptr;
  num_sb_ = 0;
  num_bb_ = 0;
}

void BitVector::bind_rank_views() {
  sb_view_ = sb_owned_.empty() ? nullptr : sb_owned_.data();
  bb_view_ = bb_owned_.empty() ? nullptr : bb_owned_.data();
}

void BitVector::push_bit(bool value) {
  if (words_view_ != nullptr) throw Error("BitVector::push_bit on a mapped vector");
  if ((size_ & 63u) == 0) {
    words_.push_back(0);
    num_words_ = words_.size();
  }
  if (value) words_.back() |= (uint64_t(1) << (63 - (size_ & 63u)));
  ++size_;
}

void BitVector::append_bits(const uint64_t* src, uint64_t nbits) {
  if (nbits == 0) return;
  if (words_view_ != nullptr) throw Error("BitVector::append_bits on a mapped vector");
  const uint64_t rem = size_ & 63u;
  const uint64_t nsrc = (nbits + 63) / 64;
  if (rem == 0) {
    words_.insert(words_.end(), src, src + nsrc);
    num_words_ = words_.size();
    size_ += nbits;
    return;
  }
  const unsigned shift = static_cast<unsigned>(64 - rem);
  words_.back() |= src[0] >> rem;
  uint64_t carry = src[0] << shift;
  const uint64_t old_words = static_cast<uint64_t>(words_.size());
  for (uint64_t i = 1; i < nsrc; ++i) {
    words_.push_back(carry | (src[i] >> rem));
    carry = src[i] << shift;
  }
  size_ += nbits;
  const uint64_t new_words = (size_ + 63) / 64;
  if (new_words > old_words + nsrc - 1) words_.push_back(carry);
  num_words_ = words_.size();
}

void BitVector::finalize() {
  num_sb_ = rank_superblock_entries(size_);
  num_bb_ = rank_basic_block_entries(size_);
  sb_owned_.assign(num_sb_, 0);
  bb_owned_.assign(num_bb_, 0);
  const uint64_t* w = words_data();
  uint64_t total = 0;
  uint64_t next_sb = 0;
  uint64_t next_bb = 0;
  for (uint64_t i = 0; i < num_words_; ++i) {
    const uint64_t pos = i << 6;
    if ((pos & (kSuperblockBits - 1)) == 0) {
      const uint64_t b = pos >> 12;
      while (next_sb <= b) sb_owned_[next_sb++] = total;
    }
    if ((pos & (kBasicBlockBits - 1)) == 0) {
      const uint64_t b = pos >> 9;
      const uint64_t base = sb_owned_[pos >> 12];
      while (next_bb <= b) bb_owned_[next_bb++] = static_cast<uint16_t>(total - base);
    }
    total += static_cast<uint64_t>(popcount64(w[i]));
  }
  while (next_sb <= num_sb_ - 1) sb_owned_[next_sb++] = total;
  while (next_bb <= num_bb_ - 1) {
    const uint64_t base = sb_owned_[(next_bb * kBasicBlockBits) >> 12];
    bb_owned_[next_bb++] = static_cast<uint16_t>(total - base);
  }
  bind_rank_views();
}

bool BitVector::bit(uint64_t i) const {
  if (i >= size_) throw QueryError("BitVector::bit: index out of range");
  return ((words_data()[i >> 6] >> (63 - (i & 63u))) & 1u) != 0;
}

uint64_t BitVector::rank1(uint64_t pos) const {
  if (sb_view_ == nullptr)
    throw Error("BitVector::rank1: rank index not built (call finalize() or load_rank())");
  if (pos > size_) throw QueryError("BitVector::rank1: position out of range");
  const uint64_t sb = pos >> 12;
  const uint64_t bb = pos >> 9;
  uint64_t r = sb_view_[sb] + bb_view_[bb];
  const unsigned bit_in_bb = static_cast<unsigned>(pos & (kBasicBlockBits - 1));
  const uint64_t base_word = bb * (kBasicBlockBits / 64);
  const uint64_t nwords = bit_in_bb >> 6;
  const uint64_t* w = words_data();
  for (uint64_t i = 0; i < nwords; ++i)
    r += static_cast<uint64_t>(popcount64(w[base_word + i]));
  const unsigned tail = bit_in_bb & 63u;
  if (tail != 0)
    r += static_cast<uint64_t>(popcount64(w[base_word + nwords] & mask_high_bits(tail)));
  return r;
}

// ---- persistence -----------------------------------------------------------

BitVector& BitVector::operator=(const BitVector& other) {
  if (this == &other) return *this;
  clear();
  size_ = other.size_;
  num_words_ = other.num_words_;
  if (num_words_ > 0) {
    const uint64_t* src = other.words_data();
    words_.assign(src, src + num_words_);
  }
  if (other.sb_view_ != nullptr) {
    sb_owned_.assign(other.sb_view_, other.sb_view_ + other.num_sb_);
    bb_owned_.assign(other.bb_view_, other.bb_view_ + other.num_bb_);
    num_sb_ = other.num_sb_;
    num_bb_ = other.num_bb_;
    bind_rank_views();
  }
  return *this;
}

void BitVector::read_file_exact(const std::string& path, std::vector<char>& out,
                                uint64_t expect) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  in.seekg(0, std::ios::end);
  const std::streamoff fsize = in.tellg();
  if (fsize < 0 || static_cast<uint64_t>(fsize) != expect)
    throw FormatError("'" + path + "': file size " + std::to_string(fsize) +
                      " does not match the expected " + std::to_string(expect) +
                      " bytes implied by the manifest");
  in.seekg(0, std::ios::beg);
  out.resize(static_cast<std::size_t>(expect));
  if (expect > 0) {
    in.read(out.data(), static_cast<std::streamsize>(expect));
    if (static_cast<uint64_t>(in.gcount()) != expect || in.bad())
      throw IoError("short read on '" + path + "'");
  }
}

FileInfo BitVector::save_bv(const std::string& path) const {
  const uint64_t total_words = bv_bytes(size_) / 8;
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path + "' for writing");
  CRC32 crc;
  const bool le = host_is_little_endian();
  const uint64_t* w = words_data();
  const uint64_t kChunkWords = 8192;  // 64 KiB buffer
  std::vector<uint8_t> buf(static_cast<std::size_t>(kChunkWords) * 8);
  for (uint64_t i = 0; i < total_words;) {
    const uint64_t chunk = std::min<uint64_t>(total_words - i, kChunkWords);
    const std::size_t bytes = static_cast<std::size_t>(chunk) * 8;
    if (le && i + chunk <= num_words_) {
      std::memcpy(buf.data(), w + i, bytes);
    } else {
      uint8_t* p = buf.data();
      for (uint64_t j = 0; j < chunk; ++j, p += 8)
        store_le64(p, (i + j) < num_words_ ? w[i + j] : 0);
    }
    out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(bytes));
    if (!out) throw OutputError("write error on '" + path + "'");
    crc.update(buf.data(), bytes);
    i += chunk;
  }
  out.flush();
  if (!out) throw OutputError("write error on '" + path + "'");
  return FileInfo{total_words * 8, crc.value()};
}

FileInfo BitVector::save_rank(const std::string& path) const {
  if (sb_view_ == nullptr)
    throw Error("BitVector::save_rank: no rank index (call finalize())");
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path + "' for writing");
  CRC32 crc;
  // Scalars are written in host byte order; the manifest records which one.
  const std::size_t sb_bytes = static_cast<std::size_t>(num_sb_) * sizeof(uint64_t);
  const std::size_t bb_bytes = static_cast<std::size_t>(num_bb_) * sizeof(uint16_t);
  std::vector<uint8_t> buf(sb_bytes + bb_bytes);
  std::memcpy(buf.data(), sb_view_, sb_bytes);
  std::memcpy(buf.data() + sb_bytes, bb_view_, bb_bytes);
  out.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()));
  out.flush();
  if (!out) throw OutputError("write error on '" + path + "'");
  crc.update(buf.data(), buf.size());
  return FileInfo{sb_bytes + bb_bytes, crc.value()};
}

void BitVector::load_bv(const std::string& path, uint64_t nbits) {
  const uint64_t expect = bv_bytes(nbits);
  std::vector<char> raw;
  read_file_exact(path, raw, expect);
  words_view_ = nullptr;
  maps_.clear();
  mapped_ = false;
  size_ = nbits;
  num_words_ = expect / 8;
  words_.assign(static_cast<std::size_t>(num_words_), 0);
  if (host_is_little_endian()) {
    if (expect > 0) std::memcpy(words_.data(), raw.data(), static_cast<std::size_t>(expect));
  } else {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(raw.data());
    for (uint64_t i = 0; i < num_words_; ++i) words_[i] = load_le64(p + 8 * i);
  }
}

void BitVector::load_rank(const std::string& path, uint64_t nbits, bool stored_little_endian) {
  const uint64_t nsb = rank_superblock_entries(nbits);
  const uint64_t nbb = rank_basic_block_entries(nbits);
  std::vector<char> raw;
  read_file_exact(path, raw, rank_bytes(nbits));
  const uint8_t* base = reinterpret_cast<const uint8_t*>(raw.data());
  const bool swap = (stored_little_endian != host_is_little_endian());
  clear_rank();
  num_sb_ = nsb;
  num_bb_ = nbb;
  sb_owned_.resize(static_cast<std::size_t>(nsb));
  bb_owned_.resize(static_cast<std::size_t>(nbb));
  for (uint64_t i = 0; i < nsb; ++i) {
    uint64_t v = 0;
    std::memcpy(&v, base + 8 * i, 8);
    sb_owned_[i] = v;
  }
  for (uint64_t j = 0; j < nbb; ++j) {
    uint16_t v = 0;
    std::memcpy(&v, base + 8 * nsb + 2 * j, 2);
    bb_owned_[j] = v;
  }
  if (swap) {
    for (uint64_t i = 0; i < nsb; ++i) {
      uint64_t v = sb_owned_[i];
      v = ((v & 0x00000000000000FFull) << 56) | ((v & 0x000000000000FF00ull) << 40) |
          ((v & 0x0000000000FF0000ull) << 24) | ((v & 0x00000000FF000000ull) << 8) |
          ((v & 0x000000FF00000000ull) >> 8) | ((v & 0x0000FF0000000000ull) >> 24) |
          ((v & 0x00FF000000000000ull) >> 40) | ((v & 0xFF00000000000000ull) >> 56);
      sb_owned_[i] = v;
    }
    for (uint64_t j = 0; j < nbb; ++j) bb_owned_[j] = bswap16(bb_owned_[j]);
  }
  bind_rank_views();
  validate_rank_index(path);
}

void BitVector::map_bv(MappedFile file, uint64_t nbits) {
  const uint64_t expect = bv_bytes(nbits);
  if (file.size() != expect)
    throw FormatError("'" + file.name() + "': mapped size " + std::to_string(file.size()) +
                      " does not match the expected " + std::to_string(expect) + " bytes");
  words_view_ = nullptr;
  words_.clear();
  size_ = nbits;
  num_words_ = expect / 8;
  if (host_is_little_endian() && file.data() != nullptr) {
    words_view_ = reinterpret_cast<const uint64_t*>(file.data());
  } else if (num_words_ > 0) {
    words_.assign(static_cast<std::size_t>(num_words_), 0);
    const uint8_t* p = file.data();
    for (uint64_t i = 0; i < num_words_; ++i) words_[i] = load_le64(p + 8 * i);
  }
  maps_.push_back(std::move(file));
  mapped_ = true;
}

void BitVector::map_rank(MappedFile file, uint64_t nbits, bool stored_little_endian) {
  const uint64_t nsb = rank_superblock_entries(nbits);
  const uint64_t nbb = rank_basic_block_entries(nbits);
  const uint64_t expect = rank_bytes(nbits);
  if (file.size() != expect)
    throw FormatError("'" + file.name() + "': mapped size " + std::to_string(file.size()) +
                      " does not match the expected " + std::to_string(expect) + " bytes");
  clear_rank();
  num_sb_ = nsb;
  num_bb_ = nbb;
  if (stored_little_endian == host_is_little_endian()) {
    sb_view_ = reinterpret_cast<const uint64_t*>(file.data());
    bb_view_ = reinterpret_cast<const uint16_t*>(file.data() + nsb * 8);
  } else {
    // Byte order differs: materialize swapped copies.
    const uint8_t* base = file.data();
    sb_owned_.resize(static_cast<std::size_t>(nsb));
    bb_owned_.resize(static_cast<std::size_t>(nbb));
    for (uint64_t i = 0; i < nsb; ++i) {
      uint64_t v = 0;
      std::memcpy(&v, base + 8 * i, 8);
      v = ((v & 0x00000000000000FFull) << 56) | ((v & 0x000000000000FF00ull) << 40) |
          ((v & 0x0000000000FF0000ull) << 24) | ((v & 0x00000000FF000000ull) << 8) |
          ((v & 0x000000FF00000000ull) >> 8) | ((v & 0x0000FF0000000000ull) >> 24) |
          ((v & 0x00FF000000000000ull) >> 40) | ((v & 0xFF00000000000000ull) >> 56);
      sb_owned_[i] = v;
    }
    for (uint64_t j = 0; j < nbb; ++j) {
      uint16_t v = 0;
      std::memcpy(&v, base + 8 * nsb + 2 * j, 2);
      bb_owned_[j] = bswap16(v);
    }
    bind_rank_views();
  }
  maps_.push_back(std::move(file));
  mapped_ = true;
  validate_rank_index(maps_.back().name());
}

void BitVector::validate_rank_index(const std::string& path) const {
  uint64_t prev = 0;
  for (uint64_t i = 0; i < num_sb_; ++i) {
    if (sb_view_[i] < prev)
      throw FormatError("'" + path + "': rank index is not monotonic");
    prev = sb_view_[i];
  }
  for (uint64_t j = 0; j < num_bb_; ++j) {
    if (bb_view_[j] > kSuperblockBits)
      throw FormatError("'" + path + "': basic-block rank entry out of range");
  }
}

}  // namespace gcrank
