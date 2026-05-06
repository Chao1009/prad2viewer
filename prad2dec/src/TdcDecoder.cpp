#include "TdcDecoder.h"

#include <cmath>
#include <limits>

using namespace tdc;

int TdcDecoder::DecodeRoc(const uint32_t *data, size_t nwords,
                          uint32_t roc_tag, TdcEventData &evt)
{
    int appended = 0;
    for (size_t i = 0; i < nwords; ++i) {
        if (evt.n_hits >= MAX_TDC_HITS) break;
        uint32_t w = data[i];
        TdcHit &h = evt.hits[evt.n_hits++];
        h.roc_tag = roc_tag;
        h.slot    = static_cast<uint8_t>((w >> 27) & 0x1F);
        h.edge    = static_cast<uint8_t>((w >> 26) & 0x1);
        h.channel = static_cast<uint8_t>((w >> 19) & 0x7F);
        h.value   = w & 0x7FFFF;
        ++appended;
    }
    return appended;
}

int TdcDecoder::DecodeReplay(const std::vector<uint32_t> &roc_tags,
                             const std::vector<uint32_t> &nwords,
                             const std::vector<uint32_t> &words,
                             TdcEventData &evt)
{
    evt.clear();
    if (roc_tags.size() != nwords.size()) return 0;

    int total = 0;
    size_t off = 0;
    for (size_t i = 0; i < roc_tags.size(); ++i) {
        size_t n = nwords[i];
        if (off + n > words.size()) break;        // truncated/corrupt entry
        total += DecodeRoc(words.data() + off, n, roc_tags[i], evt);
        off += n;
    }
    return total;
}

// --- RfTimeData --------------------------------------------------------------

float RfTimeData::nearest(const float *v, int n, float t_ref)
{
    if (n <= 0) return std::numeric_limits<float>::quiet_NaN();
    int best_i = 0;
    float best_dt = std::abs(v[0] - t_ref);
    for (int i = 1; i < n; ++i) {
        float dt = std::abs(v[i] - t_ref);
        if (dt < best_dt) { best_dt = dt; best_i = i; }
    }
    return v[best_i];
}

// --- RfTimeDecoder -----------------------------------------------------------

void RfTimeDecoder::Extract(const TdcEventData &all, RfTimeData &out,
                            uint32_t rf_roc_tag, uint8_t rf_slot,
                            uint8_t rf_ch_a,    uint8_t rf_ch_b)
{
    out.clear();
    for (int i = 0; i < all.n_hits; ++i) {
        const TdcHit &h = all.hits[i];
        if (h.roc_tag != rf_roc_tag) continue;
        if (h.slot    != rf_slot)    continue;
        if (h.edge    != 0)          continue;   // leading edges only
        if (h.channel == rf_ch_a) {
            if (out.n_a < RfTimeData::MAX_HITS_PER_CH)
                out.ns_a[out.n_a++] = static_cast<float>(h.value * TDC_LSB_NS);
        } else if (h.channel == rf_ch_b) {
            if (out.n_b < RfTimeData::MAX_HITS_PER_CH)
                out.ns_b[out.n_b++] = static_cast<float>(h.value * TDC_LSB_NS);
        }
    }
}

void RfTimeDecoder::DecodeReplay(const std::vector<uint32_t> &roc_tags,
                                 const std::vector<uint32_t> &nwords,
                                 const std::vector<uint32_t> &words,
                                 RfTimeData &out,
                                 uint32_t rf_roc_tag, uint8_t rf_slot,
                                 uint8_t rf_ch_a,    uint8_t rf_ch_b)
{
    out.clear();
    if (roc_tags.size() != nwords.size()) return;

    size_t off = 0;
    for (size_t i = 0; i < roc_tags.size(); ++i) {
        size_t n = nwords[i];
        if (off + n > words.size()) break;
        if (roc_tags[i] != rf_roc_tag) { off += n; continue; }
        for (size_t k = 0; k < n; ++k) {
            uint32_t w = words[off + k];
            uint8_t slot = static_cast<uint8_t>((w >> 27) & 0x1F);
            if (slot != rf_slot) continue;
            uint8_t edge = static_cast<uint8_t>((w >> 26) & 0x1);
            if (edge != 0) continue;
            uint8_t chan = static_cast<uint8_t>((w >> 19) & 0x7F);
            float t_ns = static_cast<float>((w & 0x7FFFF) * TDC_LSB_NS);
            if (chan == rf_ch_a) {
                if (out.n_a < RfTimeData::MAX_HITS_PER_CH)
                    out.ns_a[out.n_a++] = t_ns;
            } else if (chan == rf_ch_b) {
                if (out.n_b < RfTimeData::MAX_HITS_PER_CH)
                    out.ns_b[out.n_b++] = t_ns;
            }
        }
        off += n;
    }
}
