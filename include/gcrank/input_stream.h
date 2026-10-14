#pragma once
//
// Byte-oriented input streams (plain and gzip) and a read-only memory mapping.
//
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace gcrank {

// Abstract byte source.
class InputStream {
public:
  virtual ~InputStream();

  // Reads up to `n` bytes into `dst`; returns the number of bytes read.
  // Returns 0 only at end of input. Throws InputError on a read failure.
  virtual std::size_t read(void* dst, std::size_t n) = 0;

  // Source path, for diagnostics.
  virtual const std::string& name() const = 0;
};

// Opens `path` and picks the decompressor by sniffing the gzip magic bytes
// (0x1f 0x8b), never by looking at the file extension.
// Throws InputError when the file cannot be opened.
std::unique_ptr<InputStream> open_input(const std::string& path);

// Uncompressed file reader.
class PlainInputStream final : public InputStream {
public:
  explicit PlainInputStream(std::string path);
  ~PlainInputStream() override;

  std::size_t read(void* dst, std::size_t n) override;
  const std::string& name() const override { return name_; }

private:
  std::FILE* file_ = nullptr;
  std::string name_;
};

// gzip reader (handles concatenated gzip members).
class GzipInputStream final : public InputStream {
public:
  explicit GzipInputStream(std::string path);
  ~GzipInputStream() override;

  std::size_t read(void* dst, std::size_t n) override;
  const std::string& name() const override { return name_; }

private:
  void* gz_ = nullptr;  // gzFile, kept opaque to avoid a zlib.h include here
  std::string name_;
};

// RAII read-only memory mapping of a whole file.
//
// A default-constructed / reset MappedFile is invalid; `map()` throws IoError
// when the mapping cannot be established, so callers can fall back to reading.
// Mapping an empty file succeeds and yields data() == nullptr, size() == 0.
class MappedFile {
public:
  MappedFile() = default;
  ~MappedFile();

  MappedFile(MappedFile&& other) noexcept;
  MappedFile& operator=(MappedFile&& other) noexcept;
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  // Maps the whole file read-only. Throws IoError on failure.
  static MappedFile map(const std::string& path);

  const uint8_t* data() const { return data_; }
  std::size_t size() const { return size_; }
  // True when the object holds a usable mapping (including an empty file).
  bool valid() const { return mapped_; }
  const std::string& name() const { return name_; }

  void reset();

private:
  const uint8_t* data_ = nullptr;
  std::size_t size_ = 0;
  std::string name_;
  bool mapped_ = false;
#ifdef _WIN32
  void* file_handle_ = nullptr;  // HANDLE
  void* mapping_handle_ = nullptr;
#endif
};

}  // namespace gcrank
