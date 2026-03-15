//=============================================================================
// EvStruct                                                                  ||
// Basic information about CODA evio file format                             ||
//                                                                           ||
// Developer:                                                                ||
// Chao Peng                                                                 ||
// 09/07/2020                                                                ||
//=============================================================================
#pragma once

#include <cstdint>
#include <cstddef>


namespace evc
{

// evio bank data type
enum DataType {
    DATA_UNKNOWN32    =  (0x0),
    DATA_UINT32       =  (0x1),
    DATA_FLOAT32      =  (0x2),
    DATA_CHARSTAR8    =  (0x3),
    DATA_SHORT16      =  (0x4),
    DATA_USHORT16     =  (0x5),
    DATA_CHAR8        =  (0x6),
    DATA_UCHAR8       =  (0x7),
    DATA_DOUBLE64     =  (0x8),
    DATA_LONG64       =  (0x9),
    DATA_ULONG64      =  (0xa),
    DATA_INT32        =  (0xb),
    DATA_TAGSEGMENT   =  (0xc),
    DATA_ALSOSEGMENT  =  (0xd),
    DATA_ALSOBANK     =  (0xe),
    DATA_COMPOSITE    =  (0xf),
    DATA_BANK         =  (0x10),
    DATA_SEGMENT      =  (0x20)
};

// data word definitions
enum WordDefinition {
    BLOCK_HEADER = 0,
    BLOCK_TRAILER = 1,
    EVENT_HEADER = 2,
};

/* 32 bit bank header structure
 * ----------------------------------
 * |          length:32             |
 * ----------------------------------
 * |   tag:16   |:2| type:6 | num:8 |
 * ----------------------------------
 */
struct BankHeader
{
    uint32_t length, num, type, tag;

    BankHeader() : length(0) {}
    BankHeader(const uint32_t *buf)
    {
        length = buf[0];
        uint32_t word = buf[1];
        tag = (word >> 16) & 0xFFFF;
        type = (word >> 8) & 0x3F;
        num = (word & 0xFF);
    }

    static size_t size() { return 2; }
};

struct SegmentHeader
{
    uint32_t num, type, tag;

    SegmentHeader() {}
    SegmentHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        tag = (word >> 24) & 0xFF;
        type = (word >> 16) & 0x3F;
        num = (word & 0xFFFF);
    }

    static size_t size() { return 1; }
};

struct BlockHeader
{
    bool valid;
    uint32_t nevents, number, module, slot;

    BlockHeader() : valid(false) {}
    BlockHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == BLOCK_HEADER);
        slot = (word >> 22) & 0x1F;
        module = (word >> 18) & 0xF;
        number = (word >> 8) & 0x3FF;
        nevents = (word & 0xFF);
    }

    static size_t size() { return 1; }
};

struct BlockTrailer
{
    bool valid;
    uint32_t nwords, slot;

    BlockTrailer() : valid(false) {}
    BlockTrailer(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == BLOCK_TRAILER);
        slot = (word >> 22) & 0x1F;
        nwords = (word & 0x3FFFFF);
    }

    static size_t size() { return 1; }
};

struct EventHeader
{
    bool valid;
    uint32_t number, slot;

    EventHeader() : valid(false) {}
    EventHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == EVENT_HEADER);
        slot = (word >> 22) & 0x1F;
        number = (word & 0x3FFFFF);
    }

    static size_t size() { return 1; }
};

} // namespace evio

