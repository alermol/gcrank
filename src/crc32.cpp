#include "gcrank/crc32.h"

#include <cstdio>
#include <fstream>
#include <vector>

#include "gcrank/error.h"

namespace gcrank {

namespace {

struct Table {
  uint32_t t[256];
  Table() {
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      t[i] = c;
    }
  }
};

const Table& table() {
  static const Table kTable;  // thread-safe magic static (C++11)
  return kTable;
}

}  // namespace

void CRC32::update(const void* data, std::size_t n) {
  const uint8_t* p = static_cast<const uint8_t*>(data);
  const uint32_t* t = table().t;
  uint32_t c = value_;
  for (std::size_t i = 0; i < n; ++i) c = t[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
  value_ = c;
}

uint32_t crc32_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  CRC32 crc;
  std::vector<char> buf(1u << 18);
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got > 0) crc.update(buf.data(), static_cast<std::size_t>(got));
  }
  if (in.bad()) throw IoError("read error on '" + path + "'");
  return crc.value();
}

std::string crc32_hex(uint32_t crc) {
  char buf[9];
  std::snprintf(buf, sizeof(buf), "%08x", crc);
  return std::string(buf);
}

}  // namespace gcrank
