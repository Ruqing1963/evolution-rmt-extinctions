# evolution-rmt-extinctions

**The Rhythm of Life and Death: Random Matrix Theory Reveals Level Repulsion in Phanerozoic Mass Extinctions and Radiations**

An ecological-niche shadow imposes GOE-class repulsion on global macroevolutionary events while mixed clades randomize.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## Core Finding

The **ecological-niche-shadow hypothesis**: Earth's ecological carrying
capacity is finite. A global extinction empties niche space; a global radiation
saturates it. Either way, the biosphere needs millions of years to recharge or
reshuffle before the next comparable event — suppressing short recurrence
intervals exactly as eigenvalue repulsion suppresses small spacings. Tested
across three configurations:

| Target | Data source | n | ⟨r⟩ | CV | β | Class |
|---|---|---|---|---|---|---|
| **A: Global extinctions** | Real (Sepkoski/Bambach) | 32 | **0.618 ± 0.038** | 0.45 | 1.41 | **GOE** |
| **B: Global radiations** | Real (Sepkoski/Alroy) | 29 | **0.628 ± 0.045** | 0.50 | 1.13 | **GOE** |
| C: Mixed clades | Forward model (12 clades) | 471 | 0.359 ± 0.013 | 1.04 | 0.01 | Poisson |

Reference values: Poisson = 0.386 · GOE = 0.531 · GUE = 0.603

## The Symmetry of Creation and Destruction

The near-identical ⟨r⟩ for extinctions (0.618) and radiations (0.628) is the
central biological result: **the construction and destruction of global
biodiversity are governed by the same finite-niche repulsion.** Whether the
biosphere is emptying or filling, it cannot do so twice in immediate
succession. Both processes are limited by the same carrying capacity and the
same slow recovery dynamics.

## Repulsion Without Periodicity

Unlike the classical 26-Myr periodicity debate (Raup & Sepkoski 1984), this
analysis makes the weaker but firmer claim: extinction and radiation spacings
are **repulsive** — they follow a Wigner–Dyson distribution with strong
suppression of short intervals, rather than the exponential of a memoryless
process. This requires no clock, only a recharge time. Mean spacings (~16 Myr
extinctions, ~19 Myr radiations) match documented post-extinction recovery
intervals.

## The Biosphere as One Coupled System

The collapse to Poisson under clade mixing (Target C) is the key control. When
the coupling between lineages is removed — when life is treated as a collection
of independent clocks — repulsion vanishes. Its presence in the real global
series therefore implies that, at the level of major events, the Phanerozoic
biosphere behaves as **one coupled system sharing a finite resource.**

## Grand Synthesis — Nine Systems, One Principle

| System | Timescale | ⟨r⟩ | Study |
|---|---|---|---|
| Old Faithful geyser | ~70 min | 0.738 | Paper 5 |
| Tethyan porphyry Cu | ~1 Myr | 0.712 | Paper 4 |
| Orogenic gold | ~10 Myr | 0.678 | Paper 4 |
| Mantle plumes (×4) | ~5 Myr | 0.630 | Paper 3 |
| **Global radiations** | ~19 Myr | **0.628** | **this work** |
| **Global extinctions** | ~16 Myr | **0.618** | **this work** |
| Andean porphyry Cu | ~2 Myr | 0.601 | Paper 4 |
| Nanling W–Sn | ~5 Myr | 0.574 | Paper 4 |
| N. Atl. IRD events | ~1.5 kyr | 0.550 | Paper 5 |

**Nine physically unrelated systems — rock, fault, magma, metal, water, ice,
and the living biosphere — span only ⟨r⟩ ≈ 0.55–0.74 when a single
charge-release source is isolated, and all collapse to ⟨r⟩ ≈ 0.36–0.39 when
independent sources mix. The substrate is irrelevant; only the
charge-release-shadow structure matters.**

## Repository Structure
```
evolution-rmt-extinctions/
├── README.md · LICENSE · requirements.txt · CITATION.cff · .zenodo.json
├── paper/
│   ├── paper.tex · paper.pdf      # 12 pp.
│   └── figs/                      # PDFs embedded by LaTeX
├── code/
│   └── evolution_rmt_pipeline.py  # Three-target analysis + visualization
├── data/
│   ├── extinction_events.csv      # 32 Phanerozoic extinction peaks
│   ├── radiation_events.csv       # 29 Phanerozoic radiation peaks
│   └── mixed_clades_model.csv     # 471 forward-model events
├── figures/                       # standalone PDFs (vector)
│   ├── fig1_evolution_panel.pdf   # 1×3 spacing distributions
│   ├── fig2_phanerozoic_timeline.pdf  # event timeline
│   └── fig3_grand_synthesis.pdf   # 9-system bar chart
└── results/
    └── evolution_rmt_results.json # all statistics
```

## Method
1. **Compile** Phanerozoic extinction & radiation peaks from standard sources
2. **Local unfolding**: fit a spline to cumulative event count vs age, rescale
   spacings by local density to remove the secular (Sepkoski-curve) trend
3. **Compute** ⟨r⟩, CV, Brody β, KS tests against Poisson/GOE/GUE
4. **Control**: superpose 12 independent Poisson clades → expect Poisson

## Reproduce
```bash
pip install -r requirements.txt
cd code
python evolution_rmt_pipeline.py
```

## Six-Domain RMT Program
1. Geological boundaries (Myr) → GOE — [zenodo 20766310](https://zenodo.org/records/20766310)
2. Seismotectonic rhythms → scale-dependent — [zenodo 20768130](https://zenodo.org/records/20768130)
3. Mantle plumes (Gyr) → single-source GOE — [zenodo 20768420](https://zenodo.org/records/20768420)
4. Metallogeny (Myr) → single ore system GOE/GUE — [zenodo 20768750](https://zenodo.org/records/20768750)
5. Hydrogeology (min–kyr) → super-GUE geyser, GOE glacial floods — [github](https://github.com/Ruqing1963/hydro-rmt-geysers-floods)
6. Macroevolution (Myr) → extinction & radiation GOE, mixed-clade Poisson (**this work**)

## Honest Limitations
- Event selection involves judgment; secondary events differ between compilations (Big Five + major radiations dominate the signal)
- Phanerozoic dating uncertainties (sub-Myr to several Myr) add noise but bias *toward* Poisson, making GOE detection conservative
- Local-unfolding spline has a smoothing parameter; GOE classification stable across reasonable choices, and ⟨r⟩ is the least unfolding-sensitive diagnostic
- GOE vs GUE not resolved at n~30; robust statement is repulsion (non-Poisson), not the specific ensemble
- We document a regularity and propose a mechanism; we do not prove it against all alternatives (the Poisson collapse under mixing argues for a coupled-biosphere origin)

## Citation
```bibtex
@misc{chen2026evolution,
  author = {Chen, Ruqing},
  title  = {The Rhythm of Life and Death: Random Matrix Theory Reveals
            Level Repulsion in Phanerozoic Mass Extinctions and Radiations},
  year   = {2026},
  publisher = {GitHub},
  url    = {https://github.com/Ruqing1963/evolution-rmt-extinctions}
}
```

## License
[MIT](LICENSE). Phanerozoic event compilations derived from published
literature (Sepkoski 2002, Bambach 2006, Alroy et al. 2008, Bond & Grasby 2017).
