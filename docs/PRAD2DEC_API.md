# prad2dec — C++ API Reference

Static library for reading CODA EVIO data and decoding the PRad / PRad-II
front-end electronics. Lazy, zero-allocation in the hot path; optional
ET reader for live runs (`-DWITH_ET=ON`).

For a high-level walkthrough see [`prad2dec/README.md`](../prad2dec/README.md).
This document is the symbol reference: every public class, struct, free
function, and namespace constant exported by the library.

| Header | Namespace | Public symbols |
|---|---|---|
| [`EvChannel.h`](#evchannelh) | `evc` | `EvChannel`, `status` |
| [`EtChannel.h`](#etchannelh) | `evc` | `EtChannel` (when `WITH_ET=ON`) |
| [`EvStruct.h`](#evstructh) | `evc` | `BankHeader`, `SegmentHeader`, `TagSegmentHeader`, `EvNode`, `DataType`, `IsContainer`, `TypeName` |
| [`DaqConfig.h`](#daqconfigh) | `evc` | `DaqConfig`, `EventType`, `classify_event` |
| [`load_daq_config.h`](#load_daq_configh) | `evc` | `load_daq_config`, `load_pedestals`, `parse_hex` |
| [`Fadc250Data.h`](#fadc250datah) | `fdec` | `EventInfo`, `EventData`, `RocData`, `SlotData`, `ChannelData`, `Peak`, `Pedestal`, `DaqPeak`, `DaqWaveResult`, `PulseTemplate`, `DeconvOutput`, quality bitmasks (`Q_PEAK_*`, `Q_PED_*`, `Q_DAQ_*`, `Q_DECONV_*`) |
| [`Fadc250Decoder.h`](#fadc250decoderh) | `fdec` | `Fadc250Decoder` |
| [`Fadc250RawDecoder.h`](#fadc250rawdecoderh) | `fdec` | `Fadc250RawDecoder` |
| [`Fadc250FwAnalyzer.h`](#fadc250fwanalyzerh) | `fdec` | `Fadc250FwAnalyzer` |
| [`WaveAnalyzer.h`](#waveanalyzerh) | `fdec` | `WaveAnalyzer`, `WaveResult`, `WaveConfig` (incl. `NnlsDeconvConfig`) |
| [`PulseTemplateStore.h`](#pulsetemplatestoreh) | `fdec` | `PulseTemplateStore` |
| [`Adc1881mDecoder.h`](#adc1881mdecoderh) | `fdec` | `Adc1881mDecoder` |
| [`SspData.h`](#sspdatah) | `ssp` | `SspEventData`, `MpdData`, `ApvData`, `ApvAddress` |
| [`SspDecoder.h`](#sspdecoderh) | `ssp` | `SspDecoder` |
| [`TdcData.h`](#tdcdatah) | `tdc` | `TdcEventData`, `TdcHit` |
| [`TdcDecoder.h`](#tdcdecoderh) | `tdc` | `TdcDecoder` |
| [`VtpData.h`](#vtpdatah) | `vtp` | `VtpEventData`, `EcPeak`, `EcCluster`, `VtpBlock` |
| [`VtpDecoder.h`](#vtpdecoderh) | `vtp` | `VtpDecoder` |
| [`DscData.h`](#dscdatah) | `dsc` | `DscEventData`, `DSC2_NCH` |
| [`Dsc2Decoder.h`](#dsc2decoderh) | `dsc` | `Dsc2Decoder` |
| [`SyncData.h`](#syncdatah) | `psync` | `SyncInfo` |
| [`EpicsData.h`](#epicsdatah) | `epics` | `EpicsRecord`, `ParseEpicsText` |
| [`EpicsStore.h`](#epicsstoreh) | `epics` | `EpicsStore` |
| [`InstallPaths.h`](#installpathsh) | `prad2` | `module_dir`, `resolve_data_dir` |

---

## `EvChannel.h`

`evc::EvChannel` — the EVIO file reader. Two open modes (sequential,
random-access) with one common per-event API; lazy per-product accessors
that decode on first call and cache per sub-event.

### `evc::status`

```cpp
enum class status : int {
    failure    = -1,
    success    =  1,
    incomplete =  2,
    empty      =  3,
    eof        =  4,
};
```

### Construction and configuration

```cpp
EvChannel ch(size_t buflen = 1024 * 2000);
ch.SetConfig(const DaqConfig &cfg);   // precomputes per-product tag lists
const DaqConfig &GetConfig() const;
```

If the config's `data_banks` map is empty (legacy JSON without
`bank_structure`), default entries are synthesised from the legacy bank-tag
fields so older configs keep working.

### Open / close

```cpp
status OpenSequential(const std::string &path);
status OpenRandomAccess(const std::string &path);
status OpenAuto(const std::string &path);   // RA first, fall back to sequential
void   Close();
bool   IsRandomAccess() const;
```

### Iteration

```cpp
status Read();                              // sequential
status ReadEventByIndex(int evio_index);    // random access (0-based)
int    GetRandomAccessEventCount() const;
bool   Scan();                              // build flat node tree, classify event
EventType GetEventType() const;             // valid after Scan()
int    GetNEvents() const;                  // physics sub-events in this block
```

### Tree walk

```cpp
BankHeader GetEvHeader() const;
const std::vector<EvNode>      &GetNodes() const;
const EvNode                   &GetChild(const EvNode &n, size_t i) const;
std::vector<const EvNode*>      FindByTag(uint32_t tag) const;
const EvNode                   *FindFirstByTag(uint32_t tag) const;
const std::vector<int>         &NodesForTag(uint32_t tag) const;   // O(1)
const uint32_t *GetData (const EvNode &n) const;
const uint8_t  *GetBytes(const EvNode &n) const;
size_t          GetDataBytes(const EvNode &n) const;
const uint8_t  *GetCompositePayload(const EvNode &n, size_t &nbytes) const;
uint32_t       *GetRawBuffer();
```

### Lazy data-product accessors (preferred)

```cpp
void SelectEvent(int i) const;              // pick sub-event; clears cache if i changed
const fdec::EventInfo    &Info()  const;    // always cheap
const fdec::EventData    &Fadc()  const;    // FADC250 + ADC1881M
const ssp::SspEventData  &Gem()   const;    // SSP/MPD GEM strips
const tdc::TdcEventData  &Tdc()   const;    // V1190 timing hits
const vtp::VtpEventData  &Vtp()   const;    // VTP ECAL peaks/clusters
const psync::SyncInfo    &Sync()  const;    // run-state snapshot, persists across events
const dsc::DscEventData  &Dsc()   const;    // DSC2 scaler readout
const epics::EpicsRecord &Epics() const;    // EPICS slow-control record
```

References are invalidated by the next `Read`/`Scan`/`SelectEvent(i ≠ current)`.
`Sync()` is special: its cached `SyncInfo` persists across events and is
only refreshed when a SYNC/EPICS or control event is scanned, so physics
events have an absolute time anchor for their TI 48-bit tick delta.

### Event-anchor accessors

```cpp
int32_t  GetLastPhysicsEventNumber() const;   // -1 before first physics
uint64_t GetLastPhysicsTimestamp()   const;   // 0  before first physics
```

Slow-event consumers (EPICS, DSC2 SYNCs) use these to stamp their rows
with the most recent physics event so analysis can join trees by integer
key.

### Legacy compat path

```cpp
bool DecodeEvent(int i, fdec::EventData &evt,
                 ssp::SspEventData *ssp_evt = nullptr,
                 vtp::VtpEventData *vtp_evt = nullptr,
                 tdc::TdcEventData *tdc_evt = nullptr) const;
```

Writes into caller-owned structs without touching the lazy cache. New
code should prefer `SelectEvent` + the typed accessors.

### EPICS extraction

```cpp
std::string ExtractEpicsText() const;   // call when GetEventType() == Epics
```

Returns the raw `value  channel_name\n` text payload, or empty string.

### Debug

```cpp
void PrintTree(std::ostream &os) const;
```

`EvChannel` is non-copyable.

---

## `EtChannel.h`

`evc::EtChannel : public EvChannel` — live ET-system reader; only built
when CMake configures `-DWITH_ET=ON` (gated by `#if USE_ET`).

```cpp
EtChannel ch(size_t chunk_buf = 2000);
status Connect(const std::string &ip, int port, const std::string &et_file);
void   Disconnect();
status Open(const std::string &station);
void   Close();
status Read();
bool   IsETOpen() const;
void   AddEvFilter(std::function<bool(const BankHeader &)> &&func);
et_wrap::StationConfig &GetConfig();
```

After `Connect` + `Open`, the per-event API is identical to `EvChannel`
(`Scan` / `Info` / `Fadc` / …). Filters added via `AddEvFilter` short-
circuit `Read` for events whose top-level header doesn't pass.

---

## `EvStruct.h`

Low-level evio data structures.

`BankHeader { length, tag, type, num }`, `SegmentHeader { tag, type,
length }`, `TagSegmentHeader { tag, type, length }` — all parse a
single header from a `const uint32_t *p` argument.

`EvNode { tag, type, num, depth, parent, data_begin, data_words,
child_first, child_count }` — one node in the flat event tree built
by `Scan()`.

`enum DataType` — every evio content-type code (`DATA_BANK`,
`DATA_SEGMENT`, `DATA_TAGSEGMENT`, `DATA_COMPOSITE`, …).
`bool IsContainer(uint32_t type)` and `const char *TypeName(uint32_t)`
are inline helpers.

---

## `DaqConfig.h`

`evc::DaqConfig` — exhaustive configuration carried by `EvChannel`.
Plain struct, no JSON dependency; loaded from JSON via
[`load_daq_config.h`](#load_daq_configh).

Selected groups (every field is documented inline in the header):

- **Event-type tags** — `physics_tags`, `physics_base`, `monitoring_tags`,
  `prestart_tag`, `go_tag`, `end_tag`, `sync_tag`, `epics_tag`.
- **ADC format** — `adc_format` (`"fadc250"` or `"adc1881m"`),
  `sparsify_sigma`.
- **Bank tags** — `fadc_composite_tag`, `adc1881m_bank_tag`,
  `ti_bank_tag`, `trigger_bank_tag`, `run_info_tag`, `daq_config_tag`,
  `epics_bank_tag`, `ssp_bank_tags[]`, `fadc_raw_tag`, `tdc_bank_tag`.
- **DSC2 scaler** — `dsc_scaler.{bank_tag, slot, source ∈ {Ref,Trg,Tdc},
  channel}`; `enabled() const`.
- **FADC250 firmware emulation** — `fadc250_fw.{TET, NSB, NSA,
  MAX_PULSES, NSAT, NPED, MAXPED, CLK_NS}` (Hall-D V3 firmware names).
- **Soft analyzer** — `wave_cfg` (full `fdec::WaveConfig`).
- **TI / trigger / run-info layouts** — word-offset and shift fields for
  decoding TI bank, trigger bank (0xC000), and run-info bank (0xE10F).
- **ROC identification** — `roc_tags[]` (`{tag, name, crate, type}`),
  `ti_master_tag`.
- **SYNC / control-event format** — `sync_head_*` and `sync_control_*`
  word offsets (PRad-II defaults already populated).
- **Bank structure** — `data_banks` map: `tag → {module, product, type}`;
  used by `EvChannel`'s lazy accessors. Canonical product names exposed
  as `static constexpr` strings (`product_event_info`, `product_fadc`,
  `product_tdc`, `product_gem`, `product_vtp`, `product_sync`,
  `product_epics`, `product_daq_config`).
- **Per-channel pedestals** (ADC1881M) — `pedestals` packed by
  `pack_daq_key(crate, slot, ch)`, plus `get_pedestal(...)`.
- **Companion files** — `hycal_map_file`, `gem_map_file`,
  `pedestal_file`.

### Helpers

```cpp
bool is_physics(uint32_t tag)    const;
bool is_monitoring(uint32_t tag) const;
bool is_control(uint32_t tag)    const;
bool is_sync(uint32_t tag)       const;
bool is_epics(uint32_t tag)      const;
bool is_ssp_bank(uint32_t tag)   const;

static bool is_built_trigger_bank(uint32_t tag);   // 0xFF20-0xFF2F
static bool is_raw_trigger_bank(uint32_t tag);     // 0xFF10-0xFF1F
static bool is_trigger_bank(uint32_t tag);         // 0xFF10-0xFF4F
static bool trigger_bank_has_timestamps(uint32_t tag);
static bool trigger_bank_has_run_info(uint32_t tag);

const DataBankInfo *find_data_bank(uint32_t tag) const;
std::vector<uint32_t> banks_for_product(const std::string &product) const;
```

### `evc::EventType`

```cpp
enum class EventType : uint8_t {
    Unknown=0, Physics, Sync, Epics, Prestart, Go, End, Control,
};
EventType classify_event(uint32_t tag, const DaqConfig &cfg);
```

---

## `load_daq_config.h`

Header-only application-layer JSON loader — depends on `nlohmann/json`,
not part of the library proper.

```cpp
uint32_t evc::parse_hex(const nlohmann::json &j);   // accepts "0xFF50" or int
bool     evc::load_daq_config(const std::string &path, DaqConfig &cfg);
bool     evc::load_pedestals(const std::string &path, DaqConfig &cfg);
```

`load_daq_config` populates every field in `DaqConfig` from
`daq_config.json` including `event_tags`, `bank_tags`, `dsc_scaler`,
`fadc250_waveform.{firmware, analyzer}` (incl. `nnls_deconv`),
`ti_format`, `trigger_bank`, `run_info`, `roc_tags`, `sync_format`,
`bank_structure.data_banks`, and the companion-file pointers.

`load_pedestals` reads
`[{"crate","slot","channel","mean","rms"}, …]` into `cfg.pedestals`.

---

## `Fadc250Data.h`

Pre-allocated, flat per-event POD. No heap allocation in the event loop.
Indexed by hardware (slot, channel) for O(1) access.

### Capacity constants

| Constant | Value | |
|---|---|---|
| `MAX_SAMPLES`  | 200 | samples per channel per event |
| `MAX_CHANNELS` | 64  | channels per slot (16 FADC250, 64 ADC1881M) |
| `MAX_SLOTS`    | 32  | slot IDs 0..31 |
| `MAX_ROCS`     | 10  | ROC crates per event |
| `MAX_PEAKS`    | 8   | peaks per channel waveform |

### Event structures

`ChannelData { nsamples, samples[MAX_SAMPLES] }`.

`SlotData { present, trigger, timestamp, nchannels, channel_mask, channels[] }`.
`clear()`.

`RocData { present, tag, nslots, slots[] }`. `clear()`.

`EventInfo { type, trigger_type, trigger_bits, event_tag, event_number,
trigger_number, timestamp, run_number, unix_time }` — `trigger_type`
identifies *which* trigger fired (single value from TI event_type);
`trigger_bits` is the 32-bit FP-input mask (multiple bits possible).

`EventData { info, nrocs, roc_index, rocs[] }` with helpers `clear()` and
`const RocData *findRoc(uint32_t tag) const`.

### Soft-analyzer outputs

`Peak { height, integral, time, pos, left, right, overflow, quality }` —
`left`/`right` inclusive integration bounds; `quality` is a `Q_PEAK_*` bitmask.

`Pedestal { mean, rms, nused, quality, slope }` — `quality` is `Q_PED_*`.

### Firmware-mode outputs

`DaqPeak { pulse_id, vmin, vpeak, va, coarse, fine, time_units, time_ns,
cross_sample, peak_sample, integral, window_lo, window_hi, quality }`.

`DaqWaveResult { vnoise, npeaks, peaks[MAX_PEAKS] }`. `clear()`.

### Pulse-template / pile-up deconvolution

`PulseTemplate { tau_r_ns, tau_f_ns, is_global, grid_clk_ns,
grid[GRID_N] }` (where `GRID_OVERSAMPLE = 8`, `GRID_N = MAX_SAMPLES * 8`).

`DeconvOutput { state, n, amplitude[], height[], integral[], t0_ns[],
tau_r_ns[], tau_f_ns[], chi2_per_dof }`. `clear()`.

### Quality bitmasks

| Group | Constants |
|---|---|
| Soft peak | `Q_PEAK_GOOD = 0`, `Q_PEAK_PILED = 1<<0`, `Q_PEAK_DECONVOLVED = 1<<1` |
| Soft pedestal | `Q_PED_GOOD`, `Q_PED_NOT_CONVERGED`, `Q_PED_FLOOR_ACTIVE`, `Q_PED_TOO_FEW_SAMPLES`, `Q_PED_PULSE_IN_WINDOW`, `Q_PED_OVERFLOW`, `Q_PED_TRAILING_WINDOW` (bits 0..5) |
| Firmware peak | `Q_DAQ_GOOD`, `Q_DAQ_PEAK_AT_BOUNDARY`, `Q_DAQ_NSB_TRUNCATED`, `Q_DAQ_NSA_TRUNCATED`, `Q_DAQ_VA_OUT_OF_RANGE` (bits 0..3) |
| Deconv state (single value) | `Q_DECONV_NOT_RUN=0`, `Q_DECONV_NO_TEMPLATE=1`, `Q_DECONV_BAD_TEMPLATE=2`, `Q_DECONV_LM_NOT_CONVERGED=3`, `Q_DECONV_APPLIED=4`, `Q_DECONV_FALLBACK_GLOBAL=5` |

---

## `Fadc250Decoder.h`

```cpp
static int fdec::Fadc250Decoder::DecodeRoc(
    const uint8_t *data, size_t nbytes, RocData &roc);
```

Decode one ROC's composite FADC250 payload (format `c,i,l,N(c,Ns)`,
packed native-endian) into `roc`. Returns slot count or `-1` on error.

---

## `Fadc250RawDecoder.h`

```cpp
static int fdec::Fadc250RawDecoder::DecodeRoc(
    const uint32_t *data, size_t nwords, RocData &roc);
```

Decode the JLab FADC250 self-describing 32-bit hardware-format bank
(0xE109) when rol2 composite reformatting is skipped.

---

## `Fadc250FwAnalyzer.h`

`fdec::Fadc250FwAnalyzer` — bit-faithful firmware Mode 1/2/3 emulator.

```cpp
explicit Fadc250FwAnalyzer(const evc::DaqConfig::Fadc250FwConfig &cfg = {});
void Analyze(const uint16_t *samples, int n, float PED,
             DaqWaveResult &result) const;
evc::DaqConfig::Fadc250FwConfig cfg;
```

Reproduces TET / NSB / NSA / NSAT / NPED / MAXPED behaviour exactly as
documented in the FADC250 User's Manual. PRad-II's data has already
passed firmware TET, so the caller passes the soft analyzer's
`WaveResult.ped.mean` as `PED`.

---

## `WaveAnalyzer.h`

`fdec::WaveAnalyzer` — software waveform analysis (smoothing,
local-maxima peak finding, iterative-outlier-rejection pedestal,
integration, optional NNLS pile-up deconvolution).

### `WaveResult`

```cpp
struct WaveResult { Pedestal ped; int npeaks; Peak peaks[MAX_PEAKS]; };
```

### `WaveConfig`

Core soft-analyzer knobs (defaults shown):

| Field | Default | |
|---|---|---|
| `smooth_order`     | 2 | triangular kernel order |
| `peak_nsigma`      | 5.0 | sigma-scaled peak threshold |
| `min_peak_height`  | 10 ADC | absolute floor on detection |
| `min_peak_ratio`   | 0.3 | new peak must be ≥ this fraction of nearby peak |
| `int_tail_ratio`   | 0.1 | stop integration when below this fraction |
| `tail_break_n`     | 2 | consecutive sub-threshold samples to terminate |
| `peak_pileup_gap`  | 2 | flag `Q_PEAK_PILED` for adjacent peaks |
| `ped_nsamples`     | 30 | max pedestal-window samples |
| `ped_flatness`     | 1.0 | RMS for "flat" pedestal |
| `ped_max_iter`     | 3 | outlier-rejection iterations |
| `overflow`         | 4095 | 12-bit overflow value |
| `clk_mhz`          | 250 | clock frequency |

`NnlsDeconvConfig` (sub-struct `wave_cfg.nnls_deconv`): `enabled`,
`template_file`, `apply_to_all_peaks`, `tau_r_min/max_ns`,
`tau_f_min/max_ns`, `shape_window_factor`, `t0_window_ns`,
`amp_max_factor`, `pre_samples`, `post_samples`. See header for
detailed semantics.

### `WaveAnalyzer` methods

```cpp
explicit WaveAnalyzer(const WaveConfig &cfg = {});
void Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const;
void smooth (const uint16_t *raw, int n, float *buf) const;

void SetTemplateStore(const PulseTemplateStore *store);
const PulseTemplateStore *GetTemplateStore() const;

void SetChannelKey(int roc_tag, int slot, int channel);
void ClearChannelKey();

// Static — pulse-shape calibration utilities (used by fit_pulse_template.py)
static PulseFitResult        FitPulseShape       (const uint16_t *slice, int nslice,
                                                  int peak_idx_in_slice,
                                                  float ped, float ped_rms, float clk_ns,
                                                  float model_err_floor = 0.01f);
static PulseFitTwoTauPResult FitPulseShapeTwoTauP(/* same args */);

// Power-user / diagnostic — explicit deconv against a caller-supplied template
void Deconvolve(const uint16_t *samples, int nsamples,
                const WaveResult &wres,
                const PulseTemplate &tmpl,
                DeconvOutput &dec_out) const;

WaveConfig cfg;
```

`PulseFitResult { ok, t0_ns, tau_r_ns, tau_f_ns, peak_amp,
chi2_per_dof, n_iter }`. `PulseFitTwoTauPResult` adds a fourth `p`
parameter (rise-edge exponent).

**Thread safety:** holds per-instance mutable channel-key + template
store binding. Construct one per thread; do **not** share. The bound
`PulseTemplateStore` itself is safe to share across threads.

---

## `PulseTemplateStore.h`

`fdec::PulseTemplateStore` — owns per-type pulse templates (PbGlass,
PbWO4, LMS, Veto), loaded from the JSON written by
`analysis/pyscripts/fit_pulse_template.py`.

```cpp
bool                LoadFromFile(const std::string &path, const WaveConfig &cfg);
const PulseTemplate *Lookup(int roc_tag, int slot, int channel) const;
const PulseTemplate *type_template(const std::string &type_name) const;
bool                valid() const;
int                 n_channels_known() const;
int                 n_types_loaded()   const;
void                Clear();
```

After a successful `LoadFromFile` the store is effectively immutable;
`Lookup` is concurrent-safe. `LoadFromFile`/`Clear` are not concurrent-
safe and must be called when no `WaveAnalyzer` is using the store.

---

## `Adc1881mDecoder.h`

```cpp
struct fdec::Adc1881mDecoder {
    static constexpr uint32_t DATA_BEGIN      = 0xdc0adc00;
    static constexpr uint32_t DATA_BEGIN_MASK = 0xff0fff00;
    static constexpr uint32_t DATA_END        = 0xfabc0005;
    static constexpr uint32_t ALIGNMENT       = 0x00000000;

    static int DecodeRoc(const uint32_t *data, size_t nwords, RocData &roc);
};
```

Decodes the legacy PRad Fastbus ADC1881M format into the same `RocData`
structures FADC250 uses (single ADC value per channel as `samples[0]`,
`nsamples=1`).

---

## `SspData.h`

SSP/MPD/APV GEM readout PODs, indexed by (crate, mpd, adc_ch).

### Capacity constants

| Constant | Value |
|---|---|
| `APV_STRIP_SIZE`     | 128 |
| `SSP_TIME_SAMPLES`   | 6 |
| `MAX_APVS_PER_MPD`   | 16 |
| `MAX_MPDS`           | 64 |

### Structures

`ApvAddress { crate_id, mpd_id (fiber), adc_ch }` with `pack()` /
`operator==`.

`ApvData { addr, present, strips[128][6], nstrips, strip_mask[2], flags,
online_cm[6], has_online_cm }`. Helpers: `clear()`,
`setStrip(strip, ts, value)`, `hasStrip(strip)`. `clear()`
zero-fills `strips` (deterministic for downstream readers that don't
gate on `hasStrip`).

`MpdData { crate_id, mpd_id, present, napvs, apvs[16] }`. `clear()`.

`SspEventData { nmpds, mpds[64] }`. Helpers: `clear()`,
`findOrCreateMpd(crate, mpd)`, `findApv(crate, mpd, adc) const`.

---

## `SspDecoder.h`

```cpp
static int ssp::SspDecoder::DecodeRoc(
    const uint32_t *data, size_t nwords, int crate_id, SspEventData &evt);
```

Decodes one ROC's SSP raw bank into `evt`; `crate_id` comes from the
parent ROC bank tag mapping.

---

## `TdcData.h`

V1190 TDC PODs — flat array of hits.

| Constant | Value |
|---|---|
| `MAX_TDC_HITS` | 4096 |
| `MAX_TDC_SLOTS` | 32 |
| `MAX_TDC_CHANNELS` | 128 |

`TdcHit { roc_tag, slot, channel, edge (0=leading, 1=trailing), value }`.

`TdcEventData { n_hits, hits[MAX_TDC_HITS] }`. Helpers: `clear()`,
`countSlot(slot)`, `countChannel(slot, channel)`.

---

## `TdcDecoder.h`

```cpp
static int tdc::TdcDecoder::DecodeRoc(
    const uint32_t *data, size_t nwords, uint32_t roc_tag, TdcEventData &evt);
```

Decodes one ROC's 0xE107 payload (flat array of hits, 25 ps LSB after
rol2 normalization). Returns the number of hits appended; stops at
`MAX_TDC_HITS`.

---

## `VtpData.h`

VTP (0xE122) records — only PRad-II-relevant ECAL records are stored;
CLAS12 records are parsed past.

| Constant | Value |
|---|---|
| `MAX_EC_PEAKS`    | 512 |
| `MAX_EC_CLUSTERS` | 64 |
| `MAX_BLOCKS`      | 16 |

`EcPeak { roc_tag, inst, view, time, coord, energy }`.

`EcCluster { roc_tag, inst, time, energy, coordU, coordV, coordW }`.

`VtpBlock { roc_tag, slot, module_id, block_number, block_level, nwords,
event_number, trigger_time, has_trailer, trailer_mismatch }`.

`VtpEventData { n_peaks, n_clusters, n_blocks, peaks[], clusters[],
blocks[] }`. `clear()`.

---

## `VtpDecoder.h`

```cpp
static int vtp::VtpDecoder::DecodeRoc(
    const uint32_t *data, size_t nwords, uint32_t roc_tag, VtpEventData &evt);
```

Returns number of ECAL records (peaks + clusters) appended; `0` for
stub banks.

---

## `DscData.h`

DSC2 scaler bank decoded record.

| Constant | Value |
|---|---|
| `DSC2_NCH` | 16 |

`DscEventData { present, slot, offset, gated, ungated,
trg_gated[16], tdc_gated[16], trg_ungated[16], tdc_ungated[16],
ref_gated, ref_ungated }`. Helpers: `clear()`,
`double live_ratio() const` (returns `gated/ungated`, or `-1` when
ungated is 0).

Convention: `gated` counters are enabled while NOT busy, so live_fraction
= gated/ungated. Bank-format details are documented inline in the
header (legacy 67-word vs PRad-II rflag=1 72-word layout).

---

## `Dsc2Decoder.h`

```cpp
static bool dsc::Dsc2Decoder::DecodeBank(
    const uint32_t *data, size_t nwords,
    const evc::DaqConfig::DscScaler &cfg, DscEventData &out);

static bool dsc::Dsc2Decoder::ParsePayload(
    const uint32_t *data, size_t nwords, DscEventData &out);
```

`DecodeBank` returns `true` only when the data matches a known layout
and the slot matches `cfg.slot`; `false` otherwise (out unchanged).
`ParsePayload` is the lower-level call that fills the per-channel + ref
arrays without applying the (source, channel) selection — useful for
diagnostic tools.

---

## `SyncData.h`

`psync::SyncInfo { run_number, sync_counter, unix_time, event_tag,
run_type }` — absolute-time / run-state snapshot from SYNC and control
events. Helpers: `clear()`, `bool valid() const` (`unix_time != 0`).

The namespace is `psync`, not `sync`, to avoid collision with POSIX
`int sync(void)` from `<unistd.h>`.

---

## `EpicsData.h`

`epics::EpicsRecord { present, unix_time, sync_counter, run_number,
event_number_at_arrival, timestamp_at_arrival, channel[], value[] }` —
single-event POD. `clear()`.

```cpp
int epics::ParseEpicsText(const std::string &text, EpicsRecord &out);
```

Parses lines of `value channel_name` into the parallel `channel`/`value`
arrays. Returns the number of pairs produced. Used by both
`EvChannel::Epics()` and `EpicsStore::Feed()`.

---

## `EpicsStore.h`

`epics::EpicsStore` — run-scoped EPICS snapshot accumulator with channel
registry, value persistence (slow channels carry forward), and O(log N)
lookup by event_number.

```cpp
void Feed(int32_t event_number, uint64_t timestamp, const std::string &text);

bool GetValue(int32_t event_number, const std::string &channel, float &value) const;

struct Snapshot { int32_t event_number; uint64_t timestamp; std::vector<float> values; };
const Snapshot *FindSnapshot(int32_t event_number) const;

int                 GetChannelCount() const;
int                 GetChannelId(const std::string &name) const;
const std::string  &GetChannelName(int id) const;
const std::vector<std::string> &GetChannelNames() const;

int                GetSnapshotCount() const;
const Snapshot    &GetSnapshot(int index) const;

void Trim(int max_count);   // drop oldest snapshots
void Clear();
```

---

## `InstallPaths.h`

```cpp
std::string prad2::module_dir();
std::string prad2::resolve_data_dir(
    const char *env_name,
    std::initializer_list<const char *> rel_candidates,
    const char *compile_default);
```

Runtime data-directory resolution for shipped binaries and the `prad2py`
extension. Preference order: env var (e.g. `PRAD2_DATABASE_DIR`) → path
relative to the calling module (resolved via `dladdr` /
`GetModuleHandleExW`, so it works for both executables and the Python
extension `.so`) → build-time `DATABASE_DIR` / `RESOURCE_DIR` constant.

---

## Dependencies

- [evio](https://github.com/JeffersonLab/evio) (≥ 6.0) — required.
- [et](https://github.com/JeffersonLab/et) — optional, for `EtChannel`.

Both resolve from the Hall-B CODA installation by default; if not found,
CMake fetches from GitHub. Override with `-DEVIO_SOURCE=fetch` or
`-DET_SOURCE=fetch`.
