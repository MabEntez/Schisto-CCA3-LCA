# Data dictionary

## Field data — `data/raw/` (not distributed)

One row per child. Files: `CORNTD_recPOCCCA_2023_Tororo.csv` (n = 422),
`CORNTD_recPOCCCA_2023_Mayuge.csv` (n = 448; the export has a 449th row with no
ID, which is dropped). Both are UTF-8 with a byte-order mark.

| Column(s) | Meaning |
|---|---|
| `ID` | Participant ID |
| `Sex`, `Age`, `School`, `Date_recruitment`, `Weight`, `Height` | Demographics (not used in the analysis) |
| `Initials` | Tororo only. **Dropped on read; never carried into analysis objects.** |
| `KK_Sm1A`, `KK_Sm1B` … `KK_Sm3A`, `KK_Sm3B` | *S. mansoni* egg counts, duplicate Kato-Katz smears (A, B) on days 1–3. ×24 for eggs per gram |
| `POC_CCA_1_R1`, `_1_R2`, `_1_R3` | POC-CCA G-score (1–10), day 1, three reads (one per batch) |
| `POC_CCA_2`, `POC_CCA_3` | POC-CCA G-score, days 2 and 3 |
| `Batch_POC_CCA_urine_*` | POC-CCA batch for the matching read. Tororo: 1–3. Mayuge: lot numbers 220701075, 220902098, 221117133 (recoded to 1–3) |
| `recPOC_1_R1`, `_1_R2`, `_1_R3`, `recPOC_2`, `recPOC_3` | POC-CCA3 G-scores, as for POC-CCA |
| `Batch_recPOC_urine_*` | POC-CCA3 batch (1–3) for the matching read |
| `Box_recPOC_urine_*`, `Urine_*_stored`, `Stool_*_stored` | Logistics (not used) |
| `KK_Hkw*`, `KK_Tt*`, `KK_Asc*`, `KK_Others*` | Other helminths (not used) |
| `Pills_PZQ`, `Side_effects_treatment` | Treatment (not used) |

Missing tests are recorded as empty cells; in the models they are given a
dummy batch (4) whose precision is not interpreted.

`Calibrators_Mayuge_recPOC_CCA.csv` (field calibrators): `Batch_recPOC_CCA`,
`Box_recPOC_CCA`, `Gsc_Cal_1`, `Gsc_Cal_2`, `Gsc_Cal_3`.

## Laboratory data — `data/lab/`

`spiked_POC-CCA.csv`, `spiked_POC-CCA3.csv`: no header. Column 1 = AWA-TCA
spike (×0.03 = CCA in ng/ml); columns 2–4 = G-score (1–10) for batches 1–3.
36 rows each (12 concentrations × 3 replicate series).

## Prior — `data/priors/kints.csv`

Posterior samples of the logistic G-score curve (`k_GScore`, `intercept_Gscore`)
from the earlier POC-CCA reproducibility study; a bivariate normal fitted to
these is the prior for the curve parameters in the field models.
