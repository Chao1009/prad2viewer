# PRad-II Technical Notes

Standalone, citation-quality writeups of algorithms and data products
that ship with the PRad-II event-viewer / replay code.  Each note is
self-contained: it covers the algorithm, the C++ entry point, the
relevant config knobs, an example trace, and the figures and Python
reference implementation needed to reproduce every panel.

Each note lives in its own folder with the layout

```
<topic>/
├── <topic>.md          # the note itself
├── plots/              # PNG figures referenced by the note
└── scripts/            # NumPy + Matplotlib reference impl that
                        # regenerates the figures and prints the
                        # numeric results quoted in the note
```

## Available notes

| # | Note | Author(s) | Subject |
|---|---|---|---|
| 1 | [Software Waveform Analysis](waveform_analysis/wave_analysis.md) | Chao Peng (ANL) | `fdec::WaveAnalyzer` (median+MAD pedestal, peak finding, integration) and `fdec::Fadc250FwAnalyzer` (firmware Mode 1/2/3 emulation), worked through one example PbWO₄ pulse. |
| 2 | [HyCal Clustering](hycal_clustering/hycal_clustering.md) | Chao Peng (ANL) | `fdec::HyCalCluster` Island clustering: DFS grouping, log-weighted center-of-gravity, profile-based shower split, shower-depth correction. |
| 3 | [GEM Clustering](gem_clustering/gem_clustering.md) | Chao Peng (ANL) | `gem::GemCluster` strip-level clustering and X/Y matching: group + split + charge-weighted position, then Cartesian-with-cuts vs ADC-sorted matching. |

The author line at the top of each note is the source of truth — when
adding a new note, list the author + affiliation there as well as in
the table above.

For the firmware-mode algorithm spec (with full FADC250 manual cross-references),
see also [`docs/clas_fadc/FADC250_algorithms.md`](../clas_fadc/FADC250_algorithms.md).

## Reproducing the figures

Each note's `scripts/` folder contains a single Python file that
regenerates every figure in the corresponding `plots/` folder, using
only NumPy + Matplotlib.  From the note's directory:

```bash
cd docs/technical_notes/<topic>
python scripts/plot_<topic>.py
```

The HyCal and GEM scripts also read `database/hycal_map.json` /
`database/gem_map.json` to draw the real geometry, so they need to
be run from a checkout of this repository (the relative paths walk up
to the repo root).

## Citing a technical note

Each note is part of the PRad-II Event Viewer source tree and is
versioned with the rest of the codebase.  When citing in a paper,
talk, or memo, include the repository URL and either the commit hash
or release tag so the reader can check out the exact version you
referenced.

Use the author(s) listed in the **Author** column above (also at the
top of each note) — not the "PRad-II Collaboration" — since technical
notes are individually attributed.

**Plain-text template** — substitute `<topic>` and metadata:

> `<Author(s)>`, "`<Title>`," PRad-II Event Viewer Technical Notes,
> `<topic>`, commit `<sha>` (`<YYYY-MM-DD>`).
> https://github.com/Chao1009/prad2evviewer/blob/`<sha>`/docs/technical_notes/`<topic>`/`<topic>`.md

**BibTeX template:**

```bibtex
@misc{prad2_<topic>,
  author       = {<Last, First> and <Last, First>},
  title        = {<Title of the technical note>},
  howpublished = {PRad-II Event Viewer Technical Notes},
  year         = {<YYYY>},
  note         = {commit \texttt{<sha>}},
  url          = {https://github.com/Chao1009/prad2evviewer/blob/<sha>/docs/technical_notes/<topic>/<topic>.md},
}
```

**Worked example** (waveform analysis, replace `<sha>` / `<YYYY-MM-DD>`
with the commit you want to pin to):

```bibtex
@misc{prad2_wave_analysis,
  author       = {Peng, Chao},
  title        = {Software Waveform Analysis in {prad2dec}},
  howpublished = {PRad-II Event Viewer Technical Notes},
  year         = {2026},
  note         = {Argonne National Laboratory; commit \texttt{<sha>}},
  url          = {https://github.com/Chao1009/prad2evviewer/blob/<sha>/docs/technical_notes/waveform_analysis/wave_analysis.md},
}
```

For lab notebooks where the commit isn't important, a permalink to the
note on `main` is enough:
`https://github.com/Chao1009/prad2evviewer/blob/main/docs/technical_notes/<topic>/<topic>.md`.
