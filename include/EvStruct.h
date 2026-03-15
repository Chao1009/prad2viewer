#pragma once
//=============================================================================
// EvStruct.h — evio data structures and tree node
//=============================================================================

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace evc
{

// --- evio content type codes ------------------------------------------------
enum DataType : uint32_t {
    DATA_UNKNOWN32   = 0x0,
    DATA_UINT32      = 0x1,
    DATA_FLOAT32     = 0x2,
    DATA_CHARSTAR8   = 0x3,
    DATA_SHORT16     = 0x4,
    DATA_USHORT16    = 0x5,
    DATA_CHAR8       = 0x6,
    DATA_UCHAR8      = 0x7,
    DATA_DOUBLE64    = 0x8,
    DATA_LONG64      = 0x9,
    DATA_ULONG64     = 0xa,
    DATA_INT32       = 0xb,
    DATA_TAGSEGMENT  = 0xc,
    DATA_SEGMENT     = 0xd,
    DATA_BANK2       = 0xe,
    DATA_COMPOSITE   = 0xf,
    DATA_BANK        = 0x10,
    DATA_SEGMENT2    = 0x20,
};

inline bool IsContainer(uint32_t type)
{
    return type == DATA_BANK  || type == DATA_BANK2 ||
           type == DATA_SEGMENT || type == DATA_SEGMENT2 ||
           type == DATA_TAGSEGMENT;
}

inline const char *TypeName(uint32_t type)
{
    switch (type) {
    case DATA_UNKNOWN32:  return "UNKNOWN32";
    case DATA_UINT32:     return "UINT32";
    case DATA_FLOAT32:    return "FLOAT32";
    case DATA_CHARSTAR8:  return "STRING";
    case DATA_SHORT16:    return "SHORT16";
    case DATA_USHORT16:   return "USHORT16";
    case DATA_CHAR8:      return "CHAR8";
    case DATA_UCHAR8:     return "UCHAR8";
    case DATA_DOUBLE64:   return "DOUBLE64";
    case DATA_LONG64:     return "LONG64";
    case DATA_ULONG64:    return "ULONG64";
    case DATA_INT32:      return "INT32";
    case DATA_TAGSEGMENT: return "TAGSEG";
    case DATA_SEGMENT:    return "SEG";
    case DATA_BANK2:      return "BANK";
    case DATA_COMPOSITE:  return "COMPOSITE";
    case DATA_BANK:       return "BANK";
    case DATA_SEGMENT2:   return "SEG";
    default:              return "???";
    }
}

// --- evio header parsers (read-only, zero-copy) -----------------------------

struct BankHeader {
    uint32_t length, tag, type, num;
    BankHeader() : length(0), tag(0), type(0), num(0) {}
    BankHeader(const uint32_t *p)
        : length(p[0]),
          tag((p[1] >> 16) & 0xFFFF),
          type((p[1] >> 8) & 0x3F),
          num(p[1] & 0xFF) {}
    static constexpr size_t size() { return 2; }
    size_t data_words() const { return length >= 1 ? length - 1 : 0; }
};

struct SegmentHeader {
    uint32_t tag, type, length;   // length in words
    SegmentHeader() : tag(0), type(0), length(0) {}
    SegmentHeader(const uint32_t *p)
        : tag((p[0] >> 24) & 0xFF),
          type((p[0] >> 16) & 0x3F),
          length(p[0] & 0xFFFF) {}
    static constexpr size_t size() { return 1; }
};

struct TagSegmentHeader {
    uint32_t tag, type, length;   // length in words
    TagSegmentHeader() : tag(0), type(0), length(0) {}
    TagSegmentHeader(const uint32_t *p)
        : tag((p[0] >> 20) & 0xFFF),
          type((p[0] >> 16) & 0xF),
          length(p[0] & 0xFFFF) {}
    static constexpr size_t size() { return 1; }
};

// --- EvNode: one node in the scanned event tree -----------------------------
// Points into the raw buffer — no data is copied.
// The tree is stored as a flat vector; parent/depth give the structure.

struct EvNode {
    uint32_t tag;
    uint32_t type;
    uint32_t num;           // bank num field (0 for segments)
    int      depth;         // nesting depth (0 = top-level event)
    int      parent;        // index of parent in the node vector (-1 = root)

    // location of this node's DATA (not header) in the raw buffer
    size_t   data_begin;    // word index where data starts
    size_t   data_words;    // number of data words

    // child range in the flat node vector
    size_t   child_first;   // index of first child
    size_t   child_count;   // number of direct children
};

} // namespace evc
