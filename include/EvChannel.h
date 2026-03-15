#pragma once
//=============================================================================
// EvChannel.h — read evio events and scan the full bank tree
//=============================================================================

#include "EvStruct.h"
#include <string>
#include <vector>

namespace evc {

enum class status : int { failure = -1, success = 1, incomplete = 2, empty = 3, eof = 4 };

class EvChannel
{
public:
    EvChannel(size_t buflen = 1024 * 2000);
    virtual ~EvChannel() { Close(); }
    EvChannel(const EvChannel &) = delete;
    EvChannel &operator=(const EvChannel &) = delete;

    virtual status Open(const std::string &path);
    virtual void   Close();
    virtual status Read();

    // --- scan the current event into a flat tree of EvNodes -----------------
    // Call after Read(). Walks the entire evio structure; nothing is filtered.
    bool Scan();

    // --- accessors ----------------------------------------------------------
    BankHeader              GetEvHeader() const { return BankHeader(&buffer[0]); }
    const std::vector<EvNode> &GetNodes()  const { return nodes; }

    // children of a node
    const EvNode &GetChild(const EvNode &n, size_t i) const { return nodes[n.child_first + i]; }

    // find all nodes with a given tag (linear search, fast enough for ~100 nodes)
    std::vector<const EvNode*> FindByTag(uint32_t tag) const;

    // raw pointer into the buffer for a node's data region
    const uint32_t *GetData(const EvNode &n) const { return &buffer[n.data_begin]; }
    const uint8_t  *GetBytes(const EvNode &n) const
    {
        return reinterpret_cast<const uint8_t*>(&buffer[n.data_begin]);
    }
    size_t GetDataBytes(const EvNode &n) const { return n.data_words * sizeof(uint32_t); }

    // for composite banks: skip the tagsegment+inner-bank envelope, return payload ptr+size
    const uint8_t *GetCompositePayload(const EvNode &n, size_t &nbytes) const;

    // full raw buffer access
    uint32_t       *GetRawBuffer()       { return buffer.data(); }
    const uint32_t *GetRawBuffer() const { return buffer.data(); }

    // debug: print the scanned tree
    void PrintTree(std::ostream &os) const;

protected:
    int fHandle;
    std::vector<uint32_t> buffer;
    std::vector<EvNode>   nodes;

    // recursive scanner — returns number of words consumed
    size_t scanBank      (size_t buf_off, int depth, int parent);
    size_t scanSegment   (size_t buf_off, int depth, int parent);
    size_t scanTagSegment(size_t buf_off, int depth, int parent);
    void   scanChildren  (size_t buf_off, size_t nwords, uint32_t ptype, int depth, int parent_idx);
};

} // namespace evc
