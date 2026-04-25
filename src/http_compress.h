#pragma once
// =========================================================================
// http_compress.h — gzip helpers for the HTTP layer.
//
// The viewer's biggest endpoint (gem_apv, ~1.3 MB JSON per event) is
// highly compressible — typical 5-10× shrink.  Browsers always advertise
// Accept-Encoding: gzip; the urllib-based Python clients (gain_scanner,
// coincidence_monitor) do not, so they keep getting plain JSON.
//
// Header-only so callers don't pull in a translation unit just to gzip a
// small string; zlib does the heavy lifting.
// =========================================================================

#include <cctype>
#include <cstring>
#include <stdexcept>
#include <string>

#include <zlib.h>

namespace prad2 {

// True if the request's Accept-Encoding header lists "gzip" (case-
// insensitive token search; we don't bother with q-values).  Empty header
// or no match -> false, meaning the response goes out plain.
inline bool client_accepts_gzip(const std::string &accept_encoding)
{
    if (accept_encoding.empty()) return false;
    // Normalize to lower-case and look for a "gzip" token.  Browsers send
    // "gzip, deflate, br" or similar; the substring check is robust enough
    // because no other coding name contains "gzip".
    for (size_t i = 0; i + 4 <= accept_encoding.size(); ++i) {
        char a = std::tolower(static_cast<unsigned char>(accept_encoding[i + 0]));
        char b = std::tolower(static_cast<unsigned char>(accept_encoding[i + 1]));
        char c = std::tolower(static_cast<unsigned char>(accept_encoding[i + 2]));
        char d = std::tolower(static_cast<unsigned char>(accept_encoding[i + 3]));
        if (a == 'g' && b == 'z' && c == 'i' && d == 'p') return true;
    }
    return false;
}

// Compress `input` to gzip-wrapped deflate stream.  Returns the encoded
// bytes; throws std::runtime_error on zlib failure (only happens on OOM
// or programming bug — gzip can compress any byte string).
//
// Uses windowBits = MAX_WBITS | 16 to add the gzip wrapper that browsers
// expect (raw deflate would need windowBits negative; zlib wrapper would
// need Content-Encoding: deflate which is ambiguously implemented).
//
// Compression level 1 (Z_BEST_SPEED) is intentionally chosen — for big
// JSON payloads the speed/ratio tradeoff favours fast compression: level
// 1 typically reaches 6× ratio at ~3× the throughput of level 6.  For a
// 1.3 MB payload that is ~5 ms vs ~15 ms on a modern core.
inline std::string gzip_compress(const std::string &input, int level = 1)
{
    z_stream zs{};
    if (deflateInit2(&zs, level, Z_DEFLATED,
                     MAX_WBITS | 16,    // gzip wrapper
                     8,                  // memLevel default
                     Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("gzip_compress: deflateInit2 failed");
    }

    // Reserve a generous output buffer up front; deflateBound is the
    // worst-case size and is usually only slightly larger than input.
    std::string out;
    out.resize(deflateBound(&zs, input.size()));

    zs.next_in   = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    zs.avail_in  = static_cast<uInt>(input.size());
    zs.next_out  = reinterpret_cast<Bytef*>(out.data());
    zs.avail_out = static_cast<uInt>(out.size());

    int rc = deflate(&zs, Z_FINISH);
    deflateEnd(&zs);
    if (rc != Z_STREAM_END) {
        throw std::runtime_error("gzip_compress: deflate did not finish");
    }
    out.resize(zs.total_out);
    return out;
}

// Skip the compression cost on very small bodies — the gzip header alone
// is ~20 bytes, and HTTP framing latency dwarfs the wire savings below
// ~1 KB.  Tunable per-call by callers that know better.
inline constexpr size_t kGzipMinBytes = 1024;

} // namespace prad2
