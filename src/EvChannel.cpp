#include "EvChannel.h"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>

using namespace evc;

// --- evio C library status translation --------------------------------------
static inline status evio_status(int code)
{
    if (static_cast<unsigned>(code) == S_EVFILE_UNXPTDEOF) return status::incomplete;
    switch (code) {
    case S_SUCCESS:      return status::success;
    case EOF:            return status::eof;
    case S_EVFILE_TRUNC: return status::incomplete;
    default:             return status::failure;
    }
}

// --- construction / open / close / read -------------------------------------
EvChannel::EvChannel(size_t buflen) : fHandle(-1) { buffer.resize(buflen); }

status EvChannel::Open(const std::string &path)
{
    if (fHandle > 0) Close();
    char *cp = strdup(path.c_str()), *cm = strdup("r");
    int st = evOpen(cp, cm, &fHandle);
    free(cp); free(cm);
    return evio_status(st);
}

void EvChannel::Close() { evClose(fHandle); fHandle = -1; }

status EvChannel::Read() { return evio_status(evRead(fHandle, buffer.data(), buffer.size())); }

// === Scan: walk the entire event buffer and build the node tree =============
bool EvChannel::Scan()
{
    nodes.clear();
    BankHeader evh(&buffer[0]);
    if (evh.length + 1 > buffer.size()) return false;
    scanBank(0, 0, -1);
    return true;
}

// --- scan a BANK (2-word header) -------------------------------------------
size_t EvChannel::scanBank(size_t off, int depth, int parent)
{
    BankHeader h(&buffer[off]);
    size_t total = h.length + 1;          // header word + length words

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, h.num, depth, parent,
                     off + BankHeader::size(), h.data_words(), 0, 0});

    if (IsContainer(h.type)) {
        scanChildren(off + BankHeader::size(), h.data_words(), h.type, depth + 1, idx);
    } else if (h.type == DATA_COMPOSITE) {
        // expose the composite envelope as children (tagseg + inner bank)
        scanChildren(off + BankHeader::size(), h.data_words(), DATA_TAGSEGMENT, depth + 1, idx);
        // the second child is the inner bank — scanChildren will pick it up as a bank
        // Actually, composite = tagseg then bank, so we walk it explicitly:
        nodes[idx].child_first = 0;
        nodes[idx].child_count = 0;
        // re-scan properly:
        nodes.resize(idx + 1);  // undo any partial children

        size_t doff = off + BankHeader::size();
        size_t dwords = h.data_words();
        // tagsegment
        if (dwords >= 1) {
            size_t consumed = scanTagSegment(doff, depth + 1, idx);
            // inner bank (remainder)
            if (consumed < dwords && doff + consumed + 2 <= off + total) {
                scanBank(doff + consumed, depth + 1, idx);
            }
        }
    }
    return total;
}

// --- scan a SEGMENT (1-word header) ----------------------------------------
size_t EvChannel::scanSegment(size_t off, int depth, int parent)
{
    SegmentHeader h(&buffer[off]);
    size_t total = 1 + h.length;

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, 0, depth, parent,
                     off + 1, h.length, 0, 0});

    if (IsContainer(h.type))
        scanChildren(off + 1, h.length, h.type, depth + 1, idx);
    return total;
}

// --- scan a TAGSEGMENT (1-word header) -------------------------------------
size_t EvChannel::scanTagSegment(size_t off, int depth, int parent)
{
    TagSegmentHeader h(&buffer[off]);
    size_t total = 1 + h.length;

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, 0, depth, parent,
                     off + 1, h.length, 0, 0});

    if (IsContainer(h.type))
        scanChildren(off + 1, h.length, h.type, depth + 1, idx);
    return total;
}

// --- scan children of a container ------------------------------------------
void EvChannel::scanChildren(size_t off, size_t nwords, uint32_t ptype, int depth, int pidx)
{
    size_t first_child = nodes.size();
    size_t count = 0;
    size_t pos = 0;

    while (pos < nwords) {
        size_t consumed = 0;
        switch (ptype) {
        case DATA_BANK:
        case DATA_BANK2:
            consumed = scanBank(off + pos, depth, pidx);
            break;
        case DATA_SEGMENT:
        case DATA_SEGMENT2:
            consumed = scanSegment(off + pos, depth, pidx);
            break;
        case DATA_TAGSEGMENT:
            consumed = scanTagSegment(off + pos, depth, pidx);
            break;
        default:
            return;
        }
        if (consumed == 0) break;  // safety
        pos += consumed;
        ++count;
    }

    nodes[pidx].child_first = first_child;
    nodes[pidx].child_count = count;
}

// === accessors ==============================================================

std::vector<const EvNode*> EvChannel::FindByTag(uint32_t tag) const
{
    std::vector<const EvNode*> result;
    for (auto &n : nodes)
        if (n.tag == tag) result.push_back(&n);
    return result;
}

const uint8_t *EvChannel::GetCompositePayload(const EvNode &n, size_t &nbytes) const
{
    nbytes = 0;
    if (n.type != DATA_COMPOSITE) return nullptr;
    // composite node has 2 children: tagseg (format string) + inner bank (data)
    if (n.child_count < 2) return nullptr;
    auto &inner = nodes[n.child_first + 1];
    nbytes = inner.data_words * sizeof(uint32_t);
    return reinterpret_cast<const uint8_t*>(&buffer[inner.data_begin]);
}

// === PrintTree ==============================================================

void EvChannel::PrintTree(std::ostream &os) const
{
    for (auto &n : nodes) {
        for (int i = 0; i < n.depth; ++i) os << "  ";

        // header type label
        // determine if this was parsed as bank/seg/tagseg from depth + parent type
        os << std::setw(6) << std::left << TypeName(n.type)
           << " tag=0x" << std::hex << n.tag << std::dec << "(" << n.tag << ")"
           << " type=0x" << std::hex << n.type << std::dec
           << " num=" << n.num
           << " data=" << n.data_words << "w";

        if (n.child_count > 0)
            os << " children=" << n.child_count;

        // for leaf nodes, show a few hex words
        if (n.child_count == 0 && n.data_words > 0 && !IsContainer(n.type) && n.type != DATA_COMPOSITE) {
            os << " |";
            size_t nshow = std::min<size_t>(n.data_words, 4);
            for (size_t i = 0; i < nshow; ++i)
                os << " " << std::hex << std::setw(8) << std::setfill('0')
                   << buffer[n.data_begin + i] << std::setfill(' ') << std::dec;
            if (n.data_words > nshow) os << " ...";
        }

        // for string/char types, try to print ascii
        if ((n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) && n.data_words > 0) {
            const char *s = reinterpret_cast<const char*>(&buffer[n.data_begin]);
            size_t maxlen = n.data_words * 4;
            os << " \"";
            for (size_t i = 0; i < maxlen && s[i]; ++i) {
                if (s[i] >= 32 && s[i] < 127) os << s[i];
                else os << '.';
            }
            os << "\"";
        }

        os << "\n";
    }
}
