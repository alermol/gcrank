#include "gcrank/input_stream.h"

#include <zlib.h>

#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <utility>

#include "gcrank/error.h"

#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#else
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace gcrank {

InputStream::~InputStream() = default;

// ---- PlainInputStream ----------------------------------------------------

PlainInputStream::PlainInputStream(std::string path) : name_(std::move(path)) {
  file_ = std::fopen(name_.c_str(), "rb");
  if (file_ == nullptr)
    throw InputError("cannot open '" + name_ + "' for reading: " + std::strerror(errno));
}

PlainInputStream::~PlainInputStream() {
  if (file_ != nullptr) std::fclose(file_);
}

std::size_t PlainInputStream::read(void* dst, std::size_t n) {
  if (n == 0) return 0;
  const std::size_t got = std::fread(dst, 1, n, file_);
  if (got < n && std::ferror(file_))
    throw InputError("read error in '" + name_ + "'");
  return got;
}

// ---- GzipInputStream -----------------------------------------------------

GzipInputStream::GzipInputStream(std::string path) : name_(std::move(path)) {
  gzFile z = gzopen(name_.c_str(), "rb");
  if (z == nullptr)
    throw InputError("cannot open gzip file '" + name_ + "'");
  gz_ = z;
}

GzipInputStream::~GzipInputStream() {
  if (gz_ != nullptr) gzclose(static_cast<gzFile>(gz_));
}

std::size_t GzipInputStream::read(void* dst, std::size_t n) {
  if (n == 0) return 0;
  const std::size_t kMaxChunk = 1u << 30;  // gzread takes an unsigned length
  const unsigned want = static_cast<unsigned>(n < kMaxChunk ? n : kMaxChunk);
  const int got = gzread(static_cast<gzFile>(gz_), dst, want);
  if (got < 0) {
    int errnum = 0;
    const char* msg = gzerror(static_cast<gzFile>(gz_), &errnum);
    if (errnum == Z_ERRNO) msg = std::strerror(errno);
    throw InputError("gzip read error in '" + name_ + "': " + (msg != nullptr ? msg : "unknown"));
  }
  return static_cast<std::size_t>(got);
}

// ---- open_input ----------------------------------------------------------

std::unique_ptr<InputStream> open_input(const std::string& path) {
  std::array<unsigned char, 2> magic{};
  std::size_t got = 0;
  if (std::FILE* f = std::fopen(path.c_str(), "rb")) {
    got = std::fread(magic.data(), 1, magic.size(), f);
    std::fclose(f);
  } else {
    throw InputError("cannot open '" + path + "' for reading: " + std::strerror(errno));
  }
  if (got == magic.size() && magic[0] == 0x1f && magic[1] == 0x8b)
    return std::make_unique<GzipInputStream>(path);
  return std::make_unique<PlainInputStream>(path);
}

// ---- MappedFile ----------------------------------------------------------

MappedFile::MappedFile(MappedFile&& other) noexcept { *this = std::move(other); }

MappedFile& MappedFile::operator=(MappedFile&& other) noexcept {
  if (this != &other) {
    reset();
    data_ = std::exchange(other.data_, nullptr);
    size_ = std::exchange(other.size_, 0);
    name_ = std::move(other.name_);
    mapped_ = std::exchange(other.mapped_, false);
#ifdef _WIN32
    file_handle_ = std::exchange(other.file_handle_, nullptr);
    mapping_handle_ = std::exchange(other.mapping_handle_, nullptr);
#endif
  }
  return *this;
}

void MappedFile::reset() {
  if (!mapped_) return;
#ifdef _WIN32
  if (data_ != nullptr) ::UnmapViewOfFile(data_);
  if (mapping_handle_ != nullptr) ::CloseHandle(static_cast<HANDLE>(mapping_handle_));
  if (file_handle_ != nullptr) ::CloseHandle(static_cast<HANDLE>(file_handle_));
  file_handle_ = nullptr;
  mapping_handle_ = nullptr;
#else
  if (data_ != nullptr) ::munmap(const_cast<uint8_t*>(data_), size_);
#endif
  data_ = nullptr;
  size_ = 0;
  mapped_ = false;
}

MappedFile MappedFile::map(const std::string& path) {
  MappedFile mf;
  mf.name_ = path;
#ifndef _WIN32
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw IoError("cannot open '" + path + "' for mmap: " + std::strerror(errno));
  struct stat st;
  std::memset(&st, 0, sizeof(st));
  if (::fstat(fd, &st) != 0) {
    const int err = errno;
    ::close(fd);
    throw IoError("cannot stat '" + path + "': " + std::strerror(err));
  }
  if (st.st_size <= 0) {
    ::close(fd);
    mf.mapped_ = true;  // empty file: valid, no pages to map
    return mf;
  }
  const std::size_t length = static_cast<std::size_t>(st.st_size);
  void* p = ::mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0);
  const int err = errno;
  ::close(fd);  // the mapping stays valid after the descriptor is closed
  if (p == MAP_FAILED)
    throw IoError("mmap failed for '" + path + "': " + std::strerror(err));
  mf.data_ = static_cast<const uint8_t*>(p);
  mf.size_ = length;
  mf.mapped_ = true;
  return mf;
#else
  // Convert the UTF-8 path to UTF-16 and map read-only.
  const int wide_len = ::MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
  if (wide_len <= 0) throw IoError("cannot convert path '" + path + "' to UTF-16");
  std::wstring wide(static_cast<std::size_t>(wide_len), L'\0');
  ::MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, wide.data(), wide_len);
  HANDLE fh = ::CreateFileW(wide.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (fh == INVALID_HANDLE_VALUE) throw IoError("cannot open '" + path + "' for mmap");
  LARGE_INTEGER size;
  if (!::GetFileSizeEx(fh, &size)) {
    ::CloseHandle(fh);
    throw IoError("cannot stat '" + path + "'");
  }
  if (size.QuadPart == 0) {
    ::CloseHandle(fh);
    mf.mapped_ = true;
    return mf;
  }
  HANDLE mh = ::CreateFileMappingW(fh, nullptr, PAGE_READONLY, 0, 0, nullptr);
  if (mh == nullptr) {
    ::CloseHandle(fh);
    throw IoError("CreateFileMapping failed for '" + path + "'");
  }
  void* p = ::MapViewOfFile(mh, FILE_MAP_READ, 0, 0, 0);
  if (p == nullptr) {
    ::CloseHandle(mh);
    ::CloseHandle(fh);
    throw IoError("MapViewOfFile failed for '" + path + "'");
  }
  mf.data_ = static_cast<const uint8_t*>(p);
  mf.size_ = static_cast<std::size_t>(size.QuadPart);
  mf.file_handle_ = fh;
  mf.mapping_handle_ = mh;
  mf.mapped_ = true;
  return mf;
#endif
}

MappedFile::~MappedFile() { reset(); }

}  // namespace gcrank
