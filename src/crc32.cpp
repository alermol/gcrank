#include "gcrank/crc32.h"

#include <array>
#include <format>
#include <fstream>

#include "gcrank/error.h"

namespace gcrank {

namespace {

// CRC-32 (IEEE 802.3, reflected polynomial 0xEDB88320), computed at compile
// time so no initialization guard is needed at run time.
inline constexpr std::array<uint32_t, 256> kTable = [] {
  std::array<uint32_t, 256> table{};
  for (uint32_t i = 0; i < 256; ++i) {
    uint32_t c = i;
    for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
    table[i] = c;
  }
  return table;
}();

constexpr std::size_t kBufferBytes = 1u << 18;

}  // namespace

void CRC32::update(const void* data, std::size_t n) {
  const uint8_t* p = static_cast<const uint8_t*>(data);
  uint32_t c = value_;
  for (std::size_t i = 0; i < n; ++i) c = kTable[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
  value_ = c;
}

uint32_t crc32_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  CRC32 crc;
  std::array<char, kBufferBytes> buf;
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got > 0) crc.update(buf.data(), static_cast<std::size_t>(got));
  }
  if (in.bad()) throw IoError("read error on '" + path + "'");
  return crc.value();
}

std::string crc32_hex(uint32_t crc) {
  return std::format("{:08x}", crc);
}

}  // namespace gcrank
