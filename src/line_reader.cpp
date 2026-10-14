#include "gcrank/line_reader.h"

#include <cstring>
#include <string_view>

#include "gcrank/error.h"

namespace gcrank {

namespace {
// Length of `text` without a trailing '\r', so CRLF input reads like LF.
constexpr std::size_t without_cr(std::string_view text) {
  return !text.empty() && text.back() == '\r' ? text.size() - 1 : text.size();
}
}  // namespace

LineReader::LineReader(InputStream& in, std::size_t buffer_bytes) : in_(in) {
  if (buffer_bytes < 64) buffer_bytes = 64;
  buf_.resize(buffer_bytes);
}

bool LineReader::fill() {
  if (eof_) return false;
  begin_ = 0;
  end_ = in_.read(buf_.data(), buf_.size());
  if (end_ == 0) {
    eof_ = true;
    return false;
  }
  return true;
}

bool LineReader::next(const char** data, std::size_t* len) {
  scratch_.clear();
  for (;;) {
    if (begin_ < end_) {
      const char* base = buf_.data();
      const char* nl =
          static_cast<const char*>(std::memchr(base + begin_, '\n', end_ - begin_));
      if (nl != nullptr) {
        std::size_t line_len = static_cast<std::size_t>(nl - (base + begin_));
        const char* line = base + begin_;
        begin_ = static_cast<std::size_t>(nl - base) + 1;
        if (!scratch_.empty()) {
          // The line started in an earlier buffer; join the pieces.
          scratch_.append(line, line_len);
          line = scratch_.data();
          line_len = scratch_.size();
        }
        line_len = without_cr(std::string_view(line, line_len));
        ++line_number_;
        *data = line_len > 0 ? line : nullptr;
        *len = line_len;
        return true;
      }
      // No terminator in the rest of the buffer: keep it and read on.
      scratch_.append(base + begin_, end_ - begin_);
      begin_ = end_ = 0;
    }
    if (!fill()) {
      if (!scratch_.empty()) {
        // Unterminated last line. The pointer must stay valid until the next
        // call, so scratch_ is not cleared here.
        const std::size_t line_len = without_cr(scratch_);
        ++line_number_;
        *data = line_len > 0 ? scratch_.data() : nullptr;
        *len = line_len;
        return true;
      }
      return false;
    }
  }
}

}  // namespace gcrank
