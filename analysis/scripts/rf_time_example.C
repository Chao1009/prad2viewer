//============================================================================
// rf_time_example.C — read per-event RF time from a replayed raw file
//
// Minimal example: open a replayed *_raw.root, decode the RF reference
// per event, and print the result.
//
// The raw `events` tree carries three vector<uint32> branches per event:
//
//   tdc_roc_tags   parent ROC tag of each 0xE107 bank
//   tdc_nwords     word count per bank, parallel to tdc_roc_tags
//   tdc_words      concatenated TDC hit words
//
// (The recon tree intentionally does NOT carry these — RF-time
// reconstruction will land there as a decoded scalar in a future change.)
//
// The bit fields inside each tdc_words element (slot/edge/channel/value)
// are decoded by tdc::RfTimeDecoder::DecodeReplay() — analysis code
// should not bit-shift directly so the LSB calibration (TDC_LSB_NS) and
// PRad-II RF cabling (RF_ROC_TAG / RF_SLOT / RF_CH_A / RF_CH_B) live in
// one place (prad2dec/include/TdcData.h).
//
// In current PRad-II runs the decoder yields ~6 leading-edge hits per
// channel per event on the two RF channels (ch 0 and ch 8 of slot 16
// in ROC 0x40), spaced ~131.3 ns apart — the divided CEBAF reference.
//
// Usage (from your build dir):
//
//   root -l ../analysis/scripts/rootlogon.C
//   .x ../analysis/scripts/rf_time_example.C+( \
//       "prad_024386.00000_raw.root", 5)
//
// Args:
//   1) infile      — replayed *.root (raw or recon, auto-detected by tree)
//   2) max_events  — print at most this many events (default 5)
//
// Output: one block per event with the decoded ns arrays for both RF
// channels and three sample nearest_a() lookups so you can see the
// helper in action.
//============================================================================

#include "TdcDecoder.h"
#include "EventData.h"
#include "EventData_io.h"

#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

void rf_time_example(const char *infile, Long64_t max_events = 5)
{
    // ─── open file + auto-detect tree ─────────────────────────────────────
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[rf_time_example] cannot open " << infile << "\n";
        return;
    }
    TTree *t = dynamic_cast<TTree *>(fin->Get("events"));
    if (!t) {
        std::cerr << "[rf_time_example] no 'events' tree in " << infile
                  << " — this script reads the RAW replay (the recon tree "
                     "doesn't carry tdc_* yet)\n";
        return;
    }
    if (!t->GetBranch("tdc_words")) {
        std::cerr << "[rf_time_example] no tdc_* branches — re-replay with "
                     "the TDC-capture changes\n";
        return;
    }
    std::cout << "[rf_time_example] reading 'events' tree from " << infile << "\n";

    // ─── bind the event-id field + the three RF vector branches ───────────
    // event_num exists on both trees; the rest of the per-event payload is
    // ignored here since we only need the TDC vectors.
    int event_num = 0;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("event_num", 1);
    t->SetBranchStatus("tdc_roc_tags", 1);
    t->SetBranchStatus("tdc_nwords",   1);
    t->SetBranchStatus("tdc_words",    1);
    t->SetBranchAddress("event_num", &event_num);

    std::vector<uint32_t> *p_roc = nullptr;
    std::vector<uint32_t> *p_nw  = nullptr;
    std::vector<uint32_t> *p_w   = nullptr;
    t->SetBranchAddress("tdc_roc_tags", &p_roc);
    t->SetBranchAddress("tdc_nwords",   &p_nw);
    t->SetBranchAddress("tdc_words",    &p_w);

    // ─── per-event loop ───────────────────────────────────────────────────
    tdc::RfTimeData rf;
    Long64_t n = std::min<Long64_t>(max_events, t->GetEntries());

    std::cout << std::fixed << std::setprecision(2);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        // One call decodes the flat-triple-of-vectors representation into
        // a compact per-channel ns view (already filtered to RF_ROC_TAG /
        // RF_SLOT and leading-edge-only).
        tdc::RfTimeDecoder::DecodeReplay(*p_roc, *p_nw, *p_w, rf);

        std::cout << "\n--- event_num=" << event_num
                  << "  RF channel A (ch " << int(tdc::RF_CH_A)
                  << ", n=" << rf.n_a << "): ";
        for (int k = 0; k < rf.n_a; ++k) std::cout << " " << rf.ns_a[k];
        std::cout << " ns\n"
                  << "                       RF channel B (ch "
                  << int(tdc::RF_CH_B) << ", n=" << rf.n_b << "): ";
        for (int k = 0; k < rf.n_b; ++k) std::cout << " " << rf.ns_b[k];
        std::cout << " ns\n";

        // Demonstrate the helper that picks the nearest tick to a reference
        // time — typically the cluster time or the trigger latency you
        // want to align to.
        for (float t_ref : {200.f, 400.f, 600.f}) {
            float ta = rf.nearest_a(t_ref);
            float tb = rf.nearest_b(t_ref);
            std::cout << "    nearest at t_ref=" << t_ref << " ns: "
                      << "  A=" << ta << "  (Δ=" << (t_ref - ta) << ")"
                      << "  B=" << tb << "  (Δ=" << (t_ref - tb) << ")\n";
        }
    }

    fin->Close();
    std::cout << "\n[rf_time_example] printed " << n << " event(s).\n";
}
