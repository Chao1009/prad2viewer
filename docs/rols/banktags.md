# PRad-II EVIO Bank Tag Reference

Reference for EVIO bank tags in PRad-II data, based on run 023453 tag dump
cross-referenced with ROL sources (rol1.c, rol2.c, vtp1mpd.c) and
trigger configuration (prad_v0.trg).

## Trigger System Overview

```
FADC250 → VTP (per-crate clustering/sums) → SSP (cross-crate trigger logic)
                                                  ↓
                                            TRGBIT 0-7 via P2
                                                  ↓
                              TI/TS Front Panel inputs (bits 8-15)
                              v1495 → FP inputs (bits 24-27)
                                                  ↓
                              TS: mask + prescale + trigger table → event_type
                                                  ↓
                              Event Builder: event_tag = 0x80 + event_type
```

### FP Trigger Input Bit Map (from prad_v0.trg)

The 32-bit FP trigger inputs (TI master d[5]) have these assignments:

**SSP PRAD trigger bits (bits 8-15), routed via SSP P2 outputs:**

| Bit | Hex        | SSP P2 | Source     | Condition                              | Prescale |
|-----|------------|--------|------------|----------------------------------------|----------|
| 8   | 0x00000100 | P2OUT0 | TRGBIT 0   | Raw energy sum > 1000                  | ÷2       |
| 9   | 0x00000200 | P2OUT1 | TRGBIT 1   | >= 1 cluster, cluster E > 1000 MeV     | ÷1       |
| 10  | 0x00000400 | P2OUT2 | TRGBIT 2   | >= 2 clusters, cluster E > 1000 MeV    | ÷1       |
| 11  | 0x00000800 | P2OUT3 | TRGBIT 3   | >= 3 clusters, cluster E > 1000 MeV    | ÷1       |
| 12  | 0x00001000 | P2OUT4 | TRGBIT 4   | (disabled, mult >= 255)                | —        |
| 13  | 0x00002000 | P2OUT5 | TRGBIT 5   | (disabled, mult >= 255)                | —        |
| 14  | 0x00004000 | P2OUT6 | TRGBIT 6   | (disabled, mult >= 255)                | —        |
| 15  | 0x00008000 | P2OUT7 | TRGBIT 7   | **100 Hz pulser** (monitoring)         | ÷1       |

**v1495 external trigger bits (bits 16-27):**

| Bit | Hex        | Source   | TS prescale           |
|-----|------------|----------|-----------------------|
| 23  | 0x00800000 | v1495    | OR from SD/FADC       |
| 24  | 0x01000000 | v1495    | **LMS** — ÷1 (no prescale) |
| 25  | 0x02000000 | v1495    | alpha — ÷16385        |
| 26  | 0x04000000 | v1495    | Faraday — ÷16385      |
| 27  | 0x08000000 | v1495    | Master OR — ÷16385    |

**TS FP input mask:** `0x0F00FF00` (monitors bits 8-15 and 24-27).

### Event Tag Convention

PRad-II uses CODA single-event mode:

```
event_tag = 0x80 + TI_event_type
```

`TI_event_type` is set by `tiLoadTriggerTable(3)` based on the TS trigger
decision from the FP inputs above.

**Important:** The event_type (8-bit, in TI event header d[0] bits 31:24) and
the FP trigger inputs (32-bit, in TI master d[5]) are **different signals**:
- `event_type` = the TS trigger DECISION (which trigger type was formed)
- `d[5]` = raw FP input SNAPSHOT at trigger time (all 32 bits, including signals
  that didn't contribute to THIS trigger due to prescaling or masking)

They are related but not identical — an FP bit can be active in d[5] without
contributing to the trigger (e.g., prescaled away), and the event_type encoding
depends on the trigger table mapping.

### CODA3 Built-Trigger Tags (not used in PRad-II)

| Tag    | Description                    |
|--------|--------------------------------|
| 0xFF50 | PEB (Primary Event Builder)    |
| 0xFF58 | PEB with sync flag             |
| 0xFF70 | SEB (Secondary Event Builder)  |
| 0xFF78 | SEB with sync flag             |
| 0xFF60 | Streaming ROC Raw              |
| 0xFF62 | Streaming Data Concentrator    |
| 0xFF64 | Streaming SEB                  |
| 0xFF66 | Streaming PEB                  |

---

## Depth 0 — Top-Level Events

### Physics Events (tag = 0x80 + TI_event_type)

| Tag    | TI type | Count | Size (words)  | Has FADC | Description                        |
|--------|---------|-------|---------------|----------|------------------------------------|
| 0x00A9 | 0x29    | 796   | 240 - 71046   | Yes      | Physics trigger — FADC waveforms   |
| 0x00B0 | 0x30    | 16144 | 240 - 2105    | **No**   | Monitoring (100Hz pulser) — TI only |
| 0x00B9 | 0x39    | 773   | 24660 - 54586 | Yes      | Physics + monitoring overlap        |
| 0x00BA | 0x3A    | 2     | 240           | No       | Rare trigger                        |
| 0x00BC | 0x3C    | 443   | 240 - 305     | No       | Light trigger                       |
| 0x00FA | 0x7A    | 1     | 1063          | —        | Special / calibration               |

### Control Events (CODA2 legacy)

| Tag    | Count | Description               |
|--------|-------|---------------------------|
| 0x0011 | 2     | Prestart                  |
| 0x0012 | 2     | Go                        |
| 0x0014 | 2     | End                       |

CODA3 equivalents: 0xFFD1 (Prestart), 0xFFD2 (Go), 0xFFD4 (End).

### Other Event Types

| Tag    | Description                              |
|--------|------------------------------------------|
| 0x001F | EPICS slow control                       |
| 0x00C1 | Sync event                               |

---

## Depth 1 — ROC Banks (inside physics events)

Each physics event contains 16 children: 1 trigger bank + 1 TI master + 7 FADC crates + 7 TI slaves.

### Trigger Bank

| Tag    | Type   | Words | Description                                          |
|--------|--------|-------|------------------------------------------------------|
| 0xC000 | UINT32 | 3     | CODA trigger bank: `[event_number, event_tag, reserved]` |

### TI Master Crate

| Tag    | Type | ROC ID | Words     | Description                                |
|--------|------|--------|-----------|--------------------------------------------|
| 0x0027 | BANK | 39     | 121 - 636 | TI master — trigger info, SSP, run info    |

### HyCal FADC Crates (even tags — contain FADC waveform data)

| Tag    | Type | ROC ID | Name      | Crate | Words        |
|--------|------|--------|-----------|-------|--------------|
| 0x0080 | BANK | 128    | adchycal1 | 0     | 6 - 11768    |
| 0x0082 | BANK | 130    | adchycal2 | 1     | 6 - 12947    |
| 0x0084 | BANK | 132    | adchycal3 | 2     | 6 - 12896    |
| 0x0086 | BANK | 134    | adchycal4 | 3     | 6 - 12230    |
| 0x0088 | BANK | 136    | adchycal5 | 4     | 6 - 12845    |
| 0x008A | BANK | 138    | adchycal6 | 5     | 6 - 12896    |
| 0x008C | BANK | 140    | adchycal7 | 6     | 6 - 11871    |

Min size (6 words) = monitoring events with TI data only.
Max size (~12K words) = physics events with full FADC waveforms.

### TI Slave Crates (odd tags — TI timing data only)

| Tag    | Type | ROC ID | Name      | Crate | Words   |
|--------|------|--------|-----------|-------|---------|
| 0x0081 | BANK | 129    | tihycal1  | 0     | 6 - 89  |
| 0x0083 | BANK | 131    | tihycal2  | 1     | 6 - 90  |
| 0x0085 | BANK | 133    | tihycal3  | 2     | 6 - 90  |
| 0x0087 | BANK | 135    | tihycal4  | 3     | 6 - 90  |
| 0x0089 | BANK | 137    | tihycal5  | 4     | 6 - 90  |
| 0x008B | BANK | 139    | tihycal6  | 5     | 6 - 90  |
| 0x008D | BANK | 141    | tihycal7  | 6     | 6 - 90  |

### GEM VTP/MPD Crates (not in run 023453)

| Tag    | Type | ROC ID | Name    | Crate | ROL source  |
|--------|------|--------|---------|-------|-------------|
| 0x0031 | BANK | 49     | gemroc1 | 49    | vtp1mpd.c   |
| 0x0034 | BANK | 52     | gemroc2 | 52    | vtp1mpd.c   |

---

## Depth 2 — Data Banks (inside ROC crates)

### Banks present in run 023453

| Tag    | Type      | Count  | Words       | Parents               | ROL source                          | Description                        |
|--------|-----------|--------|-------------|-----------------------|-------------------------------------|------------------------------------|
| 0xE10A | UINT32    | 272385 | 4 - 7       | all 15 ROC crates     | `BANKOPEN(0xe10a,1,rol->pid)` → `tiReadBlock()`, stripped by rol2.c | TI hardware data         |
| 0xE101 | COMPOSITE | 11163  | 63 - 12939  | even FADC crates only | rol2.c reformats 0xE109 → 0xE101    | FADC250 composite waveforms `c,i,l,N(c,Ns)` |
| 0xE10C | UINT32    | 18159  | 102         | 0x0027 only           | `BANKOPEN(0xe10C,1,rol->pid)`       | SSP trigger processor data         |
| 0xE10F | UINT32    | 18159  | 6           | 0x0027 only           | CODA/EB                             | Run info                           |
| 0xE10E | STRING    | 16     | 81 - 6565   | 0x0011 + all ROCs     | `BANKOPEN(0xe10E,3,rol->pid)`       | DAQ config readback (first events) |

### Banks defined in ROLs but not in this run

| Tag    | Type      | ROL source                                  | Description                          |
|--------|-----------|---------------------------------------------|--------------------------------------|
| 0xE109 | UINT32    | `BANKOPEN(0xe109,1,rol->pid)` in rol1.c     | FADC250 raw hardware format (before rol2.c reformatting) |
| 0x0DEA | UINT32    | `vtpRocEbInit(VTPMPD_BANK,6,7)` in vtp1mpd.c (bank1 tag = 3562) | VTP/MPD GEM strip data (SSP bitfield format) |
| 0xE10B | UINT32    | `BANKOPEN(0xe10B,1,rol->pid)` in rol1.c     | V1190/V1290 TDC data                |
| 0xE141 | UINT32    | `BANKOPEN(0xe141,1,rol->pid)` in rol1.c     | FAV3 (FADC v3) hardware format       |
| 0xE104 | UINT32    | `BANKOPEN(0xe104,1,rol->pid)` in rol1.c     | VSCM data                           |
| 0xE105 | UINT32    | `BANKOPEN(0xe105,1,rol->pid)` in rol1.c     | DCRB / DC / Vetroc data             |
| 0xE115 | UINT32    | `BANKOPEN(0xe115,1,rol->pid)` in rol1.c     | DSC2 scaler data                    |
| 0xE112 | UINT32    | `BANKOPEN(0xe112,1,0)` in rol1.c            | HEAD bank raw format                 |
| 0xE123 | UINT32    | `BANKOPEN(0xe123,1,rol->pid)` in rol1.c     | SSP-RICH data                        |
| 0xE125 | UINT32    | `BANKOPEN(0xe125,1,ii)` in rol1.c           | Per-slot data                        |
| 0xE131 | UINT32    | `BANKOPEN(0xe131,1,rol->pid)` in rol1.c     | VFTDC data                           |
| 0xE133 | UINT32    | `BANKOPEN(0xe133,1,rol->pid)` in rol1.c     | HD (Helicity Decoder) data           |
| 0xE140 | UINT32    | `BANKOPEN(0xe140,1,0)` in rol1.c            | Special data bank                    |

---

## Depth 3 — Inside FADC Composite Banks (0xE101)

| Tag    | Type      | Count | Words       | Description                                  |
|--------|-----------|-------|-------------|----------------------------------------------|
| 0x000D | CHAR8     | 11163 | 4           | Composite format descriptor: `"c,i,l,N(c,Ns)"` |
| 0x0000 | UNKNOWN32 | 11163 | 56 - 12932  | Composite data payload (slot headers + channel samples) |

---

## Data Word Formats

### TI Data (0xE10A) — after rol2.c block header stripping

The raw TI block (from `tiReadBlock()`) includes block header/trailer.
rol2.c strips these, outputting per-event data starting with the event header.

**TI Slave (4 words):**

| Word | Content                                    | Decoder config         |
|------|--------------------------------------------|------------------------|
| d[0] | Event header: `event_type(8) \| 0x01(8) \| nwords(8)` | trigger_bits baseline: `(d[0]>>24) & 0xFF` |
| d[1] | Event number (32-bit)                      | `trigger_word = 1`     |
| d[2] | Timestamp low (32-bit)                     | `time_low_word = 2`    |
| d[3] | `evnum_high[19:16] \| ts_high[15:0]`       | `time_high_word = 3`, mask = 0xFFFF |

**TI Master (7 words) — additional words with tiSetFPInputReadout(1):**

| Word | Content                                    | Decoder config         |
|------|--------------------------------------------|------------------------|
| d[0] - d[3] | Same as slave                       |                        |
| d[4] | 8-bit trigger type pattern (TS input state) |                        |
| d[5] | **32-bit FP trigger inputs** (v1495 pattern, LMS, etc.) | `trigger_type_word = 5`, mask = 0xFFFFFFFF |
| d[6] | Additional TI data                         |                        |

The event header `event_type` (d[0] bits 31:24) and the FP trigger inputs (d[5])
are **different signals**:
- `event_type` = TS input pattern (which of 6 TS inputs fired) → determines event tag
- `d[5]` = front panel input snapshot (32 signals from v1495 etc.) → detailed trigger info

### Trigger Bank (0xC000)

| Word | Content                | Decoder config              |
|------|------------------------|-----------------------------|
| d[0] | Event number           | `trig_event_number_word = 0` |
| d[1] | Event tag (= top-level tag, e.g. 0xB0) | `trig_event_type_word = 1` |
| d[2] | Reserved (0)           |                             |

### Run Info Bank (0xE10F) — in TI master only

| Word | Content                | Decoder config              |
|------|------------------------|-----------------------------|
| d[0] | Header                 |                             |
| d[1] | Run number             | `ri_run_number_word = 1`    |
| d[2] | Event count            | `ri_event_count_word = 2`   |
| d[3] | Unix timestamp         | `ri_unix_time_word = 3`     |
| d[4] - d[5] | Additional info |                             |

### FADC250 Composite Format (0xE101)

Format string: `"c,i,l,N(c,Ns)"` — packed byte stream (little-endian).

Per slot:
- `c` (1 byte): slot ID
- `i` (4 bytes): trigger / event number
- `l` (8 bytes): 48-bit timestamp
- `N` (4 bytes): channel count

Per channel (repeated N times):
- `c` (1 byte): channel ID
- `N` (4 bytes): sample count
- `s` (2 bytes × N): ADC samples (uint16, native endian)

### FADC250 Raw Hardware Format (0xE109)

Self-describing 32-bit words with type code in bits[31:27]:

| Type (bits 31:27) | Name             | Key fields                            |
|-------------------|------------------|---------------------------------------|
| 0x00              | Block Header     | slot(5), module_id(4), block#(10), nevents(8) |
| 0x01              | Block Trailer    | slot(5), nwords(22)                   |
| 0x02              | Event Header     | slot(5), trigger#(22)                 |
| 0x03              | Trigger Time     | time(24), continuation: time_high(24) |
| 0x04              | Window Raw Data  | channel(4), width(12), then 2 ADC samples/word |
| 0x1F              | Filler           | skip                                  |

Sample data words (continuation after type 0x04 header):
- Bits 28:16 = ADC sample i (13 bits)
- Bits 12:0 = ADC sample i+1 (13 bits)

### SSP/MPD Bitfield Format (0xE10C, 0x0DEA)

Self-describing 32-bit words. Bit 31 = type-defining flag, bits 30:27 = type tag.

| Type | Name              | Key fields                              |
|------|-------------------|-----------------------------------------|
| 0    | Block Header      | slot(5), module_id(4), block#(10), nevents(8) |
| 1    | Block Trailer     | words_in_block(22), slot(5)             |
| 2    | Event Header      | trigger_number(27)                      |
| 3    | Trigger Time      | low(24) + continuation high(24)         |
| 5    | MPD Frame / APV   | Defining: mpd_id(5), fiber(6), flags(5). Continuation: 3 words per strip with 6 time samples (13-bit signed each) |
| 4, 6-12 | MPD Timestamp  | timestamp_fine(8), coarse(16+24), event_count(20) |
| 0xD  | MPD Debug         | Online common mode: 6 × 13-bit CM values across 3 words |
| 0xE  | Data Not Valid     | skip                                    |
| 0xF  | Filler            | skip                                    |

---

## Event Structure Diagram

```
[Physics Event: tag = 0x80 + TI_event_type] (depth 0)
 ├── [0xC000 UINT32, 3w] — trigger bank
 ├── [0x0027 BANK] — TI master crate
 │    ├── [0xE10A UINT32, 7w] — TI data (trigger#, timestamp, FP trigger bits)
 │    ├── [0xE10C UINT32, 102w] — SSP trigger processor data
 │    ├── [0xE10E STRING] — DAQ config string (first events only)
 │    └── [0xE10F UINT32, 6w] — run info (run#, unix time)
 ├── [0x0080 BANK] — HyCal FADC crate 1
 │    ├── [0xE10A UINT32, 4w] — TI data
 │    └── [0xE101 COMPOSITE] — FADC250 waveforms (physics triggers only)
 ├── [0x0081 BANK] — TI slave for crate 1
 │    └── [0xE10A UINT32, 4w] — TI data only
 ├── ... (7 FADC crates 0x80-0x8C + 7 TI slaves 0x81-0x8D)
 └── [0x008D BANK] — TI slave for crate 7
      └── [0xE10A UINT32, 4w]

[Control Event: tag = 0x0011/0x0012/0x0014] (depth 0)
 └── [0xE10E STRING] — DAQ config / run parameters
```
