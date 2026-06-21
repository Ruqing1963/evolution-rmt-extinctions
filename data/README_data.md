# Data

## extinction_events.csv
32 Phanerozoic global extinction peaks spanning 510–5 Ma, including the Big Five.
Compiled from Sepkoski (2002), Bambach (2006), Raup & Sepkoski (1982), and
Bond & Grasby (2017), with ages updated to the ICS 2023 timescale.

**Columns:**
- `event_number`: sequential index (old to young)
- `age_ma`: age in millions of years before present
- `event_name`: stratigraphic/event name
- `big_five`: yes/no flag for the five largest mass extinctions

The Big Five: End-Ordovician (444 Ma), Late Devonian/Kellwasser (372 Ma),
End-Permian (252 Ma), End-Triassic (201 Ma), End-Cretaceous (66 Ma).

## radiation_events.csv
29 Phanerozoic global radiation (origination/diversification) peaks spanning
530–5 Ma, including the Cambrian Explosion and the Great Ordovician
Biodiversification Event (GOBE). Compiled from Sepkoski (2002), Alroy et al.
(2008), and Servais et al. (2009).

**Columns:**
- `event_number`: sequential index (old to young)
- `age_ma`: age in millions of years before present
- `event_name`: radiation event name

## mixed_clades_model.csv
Forward-model output: 471 background-turnover events from 12 ecologically
isolated clades, each an independent Poisson process (rates 0.04–0.12
events/Myr) over the Phanerozoic 540-Myr window. Seed = 1859 (Darwin's
*Origin of Species*). Used to validate the spectral superposition theorem.

**Columns:**
- `event_number`: sequential index
- `age_ma`: event age in Ma (merged and sorted across all 12 clades)

## Note on unfolding
Phanerozoic event density is non-stationary (the Sepkoski curve). The analysis
applies *local unfolding* — a spline fit to cumulative event count vs age,
rescaling spacings by local density — before computing spacing statistics. See
`code/evolution_rmt_pipeline.py` for the implementation.
