#pragma once
//
// Buffered line reader over an InputStream.
//
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "gcrank/input_stream.h"

namespace gcrank {

// Splits an input stream into lines.
//
// A line ends at '\n'; a single '\r' immediately before it is stripped, so
// both LF and CRLF files are read correctly. Lines may span the internal
// buffer; the reported pointer stays valid until the next call to next().
class LineReader {
public:
  explicit LineReader(InputStream& in, std::size_t buffer_bytes = 1u << 16);

  // Advances to the next line. Returns false at end of input.
  // On success, `*data` points to the line content (without the terminator)
  // and `*len` is its length. `*data` may be nullptr when the line is empty.
  bool next(const char** data, std::size_t* len);

  // 1-based number of the line last returned.
  uint64_t line_number() const { return line_number_; }

private:
  bool fill();

  InputStream& in_;
  std::vector<char> buf_;
  std::size_t begin_ = 0;
  std::size_t end_ = 0;
  bool eof_ = false;
  std::string scratch_;  // only used when a line spans buffer boundaries
  uint64_t line_number_ = 0;
};

}  // namespace gcrank
