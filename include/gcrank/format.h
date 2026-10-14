#pragma once
//
// On-disk format contract for `.gcrank` indexes, format_version "1.0.0".
//
// Layout of `<NAME>.gcrank/`:
//
//   manifest.json     canonical JSON (sorted keys, 2-space indent, LF)
//   contigs.tsv       <name>\t<length>\t<offset>  (LF, no header, 0-based)
//   genome.bv         GC bitvector (1 = unambiguous G/C)
//   genome.rb         rank index for genome.bv
//   valid.bv          valid-base bitvector (1 = unambiguous A/C/G/T)
//   valid.rb          rank index for valid.bv
//   <file>.crc32      sidecar checksums, "<8 lowercase hex digits>\n"
//
// Bit order (bit_order = "msb0"): bit `i` of a bitvector is bit `(7 - i%8)` of
// byte `i/8` of the file. The stream of bytes is grouped into 64-bit words for
// the in-memory representation with the first byte in the least significant
// position; this grouping is fixed by the format and is host-independent.
//
// genome.bv / valid.bv: `bv_bytes(nbits)` bytes, zero padded to a whole number
// of superblocks (4096 bits). The logical size (total_bases) lives in the
// manifest. Ambiguous input (N, U, R, Y, S, W, K, M, B, D, H, V) sets neither
// bit: those positions are neither GC nor valid.
//
// genome.rb / valid.rb: multilbyte scalars stored in the byte order recorded
// by the manifest's `endianness` field; readers byte-swap when the host order
// differs. The layout is
//   (num_superblocks(nbits) + 1) x uint64   absolute rank at i*4096
//   (num_basic_blocks(nbits) + 1) x uint16  rank within the superblock at j*512
// The extra trailing entry of each array is the rank at the end of the vector,
// so rank1(pos) is defined for every pos in [0, nbits].
//
#include <cstdint>
#include <string>

#include "gcrank/error.h"

namespace gcrank {

// ---- versions ------------------------------------------------------------

// Version of the on-disk format (the public contract).
inline constexpr const char* kFormatVersion = "1.0.0";
// Version of the writer (library/builder); changes freely across releases.
inline constexpr const char* kBuilderVersion = "1.0.0";

// ---- index layout --------------------------------------------------------

inline constexpr uint32_t kSuperblockBits = 4096;
inline constexpr uint32_t kBasicBlockBits = 512;
inline constexpr const char* kBitOrderMsb0 = "msb0";

inline constexpr const char* kManifestFileName = "manifest.json";
inline constexpr const char* kContigsFileName = "contigs.tsv";
inline constexpr const char* kGenomeBvFileName = "genome.bv";
inline constexpr const char* kGenomeRbFileName = "genome.rb";
inline constexpr const char* kValidBvFileName = "valid.bv";
inline constexpr const char* kValidRbFileName = "valid.rb";
inline constexpr const char* kCrcSuffix = ".crc32";

inline constexpr const char* kEndiannessLittle = "little";
inline constexpr const char* kEndiannessBig = "big";

// Data files covered by checksums (manifest `files` and .crc32 sidecars).
inline constexpr const char* kDataFiles[] = {
    kContigsFileName, kGenomeBvFileName, kGenomeRbFileName,
    kValidBvFileName, kValidRbFileName,
};

// ---- sizes ---------------------------------------------------------------

inline constexpr uint64_t num_superblocks(uint64_t nbits) {
  return (nbits + kSuperblockBits - 1) / kSuperblockBits;
}
inline constexpr uint64_t num_basic_blocks(uint64_t nbits) {
  return (nbits + kBasicBlockBits - 1) / kBasicBlockBits;
}
// Physical size of a .bv file for `nbits` logical bits (zero padded).
inline constexpr uint64_t bv_bytes(uint64_t nbits) {
  return num_superblocks(nbits) * (kSuperblockBits / 8);
}
inline constexpr uint64_t rank_superblock_entries(uint64_t nbits) {
  return num_superblocks(nbits) + 1;
}
inline constexpr uint64_t rank_basic_block_entries(uint64_t nbits) {
  return num_basic_blocks(nbits) + 1;
}
// Physical size of a .rb file for `nbits` logical bits.
inline constexpr uint64_t rank_bytes(uint64_t nbits) {
  return rank_superblock_entries(nbits) * sizeof(uint64_t) +
         rank_basic_block_entries(nbits) * sizeof(uint16_t);
}

// Size and checksum of one data file (manifest `files` entry).
struct FileInfo {
  uint64_t size = 0;
  uint32_t crc32 = 0;
};

// ---- host properties -----------------------------------------------------

// True when the host stores multi-byte integers least significant byte first.
inline bool host_is_little_endian() {
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
  return __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__;
#elif defined(_WIN32)
  return true;  // every Windows target supported by gcrank is little-endian
#else
  const uint16_t probe = 1;
  return *reinterpret_cast<const uint8_t*>(&probe) == 1;
#endif
}

inline const char* host_endianness_name() {
  return host_is_little_endian() ? kEndiannessLittle : kEndiannessBig;
}

// ---- base semantics ------------------------------------------------------
//
// Semantics of `genome.bv` and `valid.bv`, and of the FASTA parser:
//   genome.bv bit = 1 for an unambiguous G/C
//   valid.bv  bit = 1 for an unambiguous A/C/G/T (U counts as T)
//
// Accepted sequence bytes (upper and lower case are equivalent; the parser
// normalizes to upper case):
//
//   A C G T   unambiguous DNA. A/C/G/T sets the valid bit; G/C also sets the
//             GC bit.
//   U         uracil, the RNA spelling of T. It is normalized to T and is
//             therefore valid but never GC (this is what `bedtools nuc` does
//             as well).
//   N         unknown base.
//   R Y S W K M B D H V
//             IUPAC ambiguity codes (R = A/G, Y = C/T, S = G/C, W = A/T,
//             K = G/T, M = A/C, B = C/G/T, D = A/G/T, H = A/C/T, V = A/C/G).
//
// Ambiguous input sets neither bit: it is neither GC nor valid, so it never
// contributes to gc_count and gc_fraction(remove_ambiguous = true) removes it
// from the denominator. Note that a code which merely *may* be G/C (S, B, V)
// is not counted as GC, so gc_count is a lower bound for windows that contain
// ambiguity codes.
//
// Whitespace inside sequence lines is ignored; any other byte (including gaps
// such as '-', '.' and '*') is an InputError.

enum BaseClass : uint8_t {
  kBaseInvalid = 0,
  kBaseA = 1,
  kBaseC = 2,
  kBaseG = 3,
  kBaseT = 4,
  kBaseAmbiguous = 5,  // N and the IUPAC ambiguity codes
  kBaseWhitespace = 6,
};

// The IUPAC ambiguity codes and N; upper case (the lower case letters are
// mapped to the same class).
inline constexpr const char* kAmbiguousBasesUpper = "NRYSWKMBDHV";
inline constexpr const char* kAmbiguousBasesLower = "nryswkmbdhv";

struct BaseTable {
  uint8_t v[256] = {};
  constexpr BaseTable() {
    v[static_cast<unsigned>('A')] = kBaseA;
    v[static_cast<unsigned>('C')] = kBaseC;
    v[static_cast<unsigned>('G')] = kBaseG;
    v[static_cast<unsigned>('T')] = kBaseT;
    v[static_cast<unsigned>('a')] = kBaseA;
    v[static_cast<unsigned>('c')] = kBaseC;
    v[static_cast<unsigned>('g')] = kBaseG;
    v[static_cast<unsigned>('t')] = kBaseT;
    // Uracil is the RNA spelling of thymine.
    v[static_cast<unsigned>('U')] = kBaseT;
    v[static_cast<unsigned>('u')] = kBaseT;
    for (int i = 0; kAmbiguousBasesUpper[i] != '\0'; ++i)
      v[static_cast<unsigned char>(kAmbiguousBasesUpper[i])] = kBaseAmbiguous;
    for (int i = 0; kAmbiguousBasesLower[i] != '\0'; ++i)
      v[static_cast<unsigned char>(kAmbiguousBasesLower[i])] = kBaseAmbiguous;
    v[static_cast<unsigned>(' ')] = kBaseWhitespace;
    v[static_cast<unsigned>('\t')] = kBaseWhitespace;
    v[static_cast<unsigned>('\v')] = kBaseWhitespace;
    v[static_cast<unsigned>('\f')] = kBaseWhitespace;
    v[static_cast<unsigned>('\r')] = kBaseWhitespace;
  }
};

inline constexpr BaseTable kBaseTable{};

// Upper-case normalization of accepted sequence bytes, with U (uracil) folded
// onto T. Applied to every base the parser hands out.
struct UpperTable {
  uint8_t v[256] = {};
  constexpr UpperTable() {
    for (int c = 0; c < 256; ++c) v[c] = static_cast<uint8_t>(c);
    for (int c = 'a'; c <= 'z'; ++c) v[c] = static_cast<uint8_t>(c - 32);
    v[static_cast<unsigned>('U')] = static_cast<uint8_t>('T');
    v[static_cast<unsigned>('u')] = static_cast<uint8_t>('T');
  }
};

inline constexpr UpperTable kUpperTable{};

// Classification lookup for one byte of sequence data.
inline constexpr uint8_t base_class(char c) {
  return kBaseTable.v[static_cast<unsigned char>(c)];
}
// Upper-case (and U -> T) normalization for the parsed sequence.
inline constexpr char normalized_base(char c) {
  return static_cast<char>(kUpperTable.v[static_cast<unsigned char>(c)]);
}
// True for an unambiguous G/C (the genome.bv bit).
inline constexpr bool base_is_gc(uint8_t cls) {
  return cls == kBaseG || cls == kBaseC;
}
// True for an unambiguous A/C/G/T/U (the valid.bv bit).
inline constexpr bool base_is_valid(uint8_t cls) {
  return cls >= kBaseA && cls <= kBaseT;
}
// True for N and the IUPAC ambiguity codes.
inline constexpr bool base_is_ambiguous(uint8_t cls) { return cls == kBaseAmbiguous; }

// ---- version parsing -----------------------------------------------------

// Returns the major version of a "MAJOR.MINOR.PATCH" string.
// Throws FormatError when the string is not a version number.
inline int parse_major_version(const std::string& version) {
  std::size_t i = 0;
  if (i >= version.size() || version[i] < '0' || version[i] > '9')
    throw FormatError("malformed version string '" + version + "'");
  long major = 0;
  while (i < version.size() && version[i] >= '0' && version[i] <= '9') {
    major = major * 10 + (version[i] - '0');
    if (major > 1000000) throw FormatError("malformed version string '" + version + "'");
    ++i;
  }
  if (i >= version.size() || version[i] != '.')
    throw FormatError("malformed version string '" + version + "'");
  return static_cast<int>(major);
}

}  // namespace gcrank
