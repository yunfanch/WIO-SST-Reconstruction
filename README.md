# WIO-SST-Reconstruction

Multi-method reconstruction of Western Indian Ocean (WIO) sea surface temperature (SST) using coral geochemical proxies.

This repository contains the analysis code for:

> Chen, Y., Cole, J., Stevenson, S., Feng, Y., Barnett, H., Tanuputri, J., Dyez, K., McClanahan, T. *Indo-Pacific Interactions in the Industrial Era: Insights from Coral Records*. (manuscript in review)

## Overview

We combine four new coral δ¹⁸O chronologies from Tanzania with seven published coral records from the WIO to reconstruct 160 years of regional SST (1846–2006). Four independent reconstruction methods are applied — Composite Plus Scale (CPS), Principal Component Analysis (PCA), Principal Component Regression (PCR), and Regularized Expectation Maximization (RegEM) — and validated against three instrumental SST products (HadISST, ERSST, COBE-SST2).

## Repository Structure

```
WIO-SST-Reconstruction/
├── Reconstruction/         # Core reconstruction methods
│   ├── cps_weighted.m
│   ├── pca_Reconstruction.m
│   ├── pcr_Reconstruction.m
│   ├── pcrsse.m
│   ├── regem_perform.m
│   ├── regem_Reconstruction.m
│   ├── nested_Reconstruction.m
│   └── performReconstruction.m
├── Validation/              # Skill metrics
│   ├── RE.m
│   └── CE.m
├── Utils/                   # Helper / pre-processing functions
│   ├── aggregate.m
│   └── scale.m
├── .gitignore
├── LICENSE
└── README.md
```

## Dependencies

- MATLAB (developed/tested with R2023b or later; no toolboxes beyond base MATLAB and the Statistics and Machine Learning Toolbox are required)
- **[RegEM](https://github.com/tapios/RegEM)** (Schneider, 2001) — required for `regem_perform.m` / `regem_Reconstruction.m`. This package is **not redistributed here** because it is released under the GPL-3.0 license. Please clone it separately and add it to your MATLAB path:

  ```bash
  git clone https://github.com/tapios/RegEM.git
  ```

  ```matlab
  addpath('/path/to/RegEM')
  ```

## Reconstruction Methods

| Method | Script | Description |
| --- | --- | --- |
| CPS | `cps_weighted.m` | Composite-plus-scale; proxies weighted by correlation with target SST |
| PCA | `pca_Reconstruction.m` | First principal component of the proxy network, sign-adjusted to target |
| PCR | `pcr_Reconstruction.m` | Linear regression of target SST on retained PCs (selected via 10-fold CV) |
| RegEM | `regem_Reconstruction.m`, `regem_perform.m` | Regularized EM with truncated total least squares (TTLS), via Schneider (2001) |
| — | `nested_Reconstruction.m` | Nested approach extending the reconstruction back in time as shorter proxy series drop out |
| — | `performReconstruction.m` | Top-level driver that runs all four methods over the nested proxy network |

Skill is assessed using Pearson correlation (r), reduction of error (RE), and coefficient of efficiency (CE), computed in `Validation/RE.m` and `Validation/CE.m` over independent calibration and verification windows.

## Data Availability

Raw coral and instrumental data are not included in this repository. They are available from:

- New Tanzania coral δ¹⁸O records: [NOAA NCEI Paleoclimatology Database](https://www.ncei.noaa.gov/access/paleo-search/study/44219)
- Published coral records: [NOAA NCEI Paleoclimatology Database](https://www.ncei.noaa.gov/products/paleoclimatology)
- HadISSTv1.1: [Met Office Hadley Centre](https://www.metoffice.gov.uk/hadobs/hadisst/)
- ERSSTv5 / OISSTv2.1: [NOAA NCEI](https://www.ncei.noaa.gov/)
- COBE-SST2, ICOADS, 20CRv3: [NOAA PSL](https://psl.noaa.gov/)

## Usage

1. Clone this repository and add all subfolders to your MATLAB path.
2. Clone and add [RegEM](https://github.com/tapios/RegEM) to your path if running the RegEM reconstruction.
3. Download the proxy and instrumental SST data listed above.
4. Run `performReconstruction.m`, adjusting input paths to point to your local data.

## Citation

If you use this code, please cite the manuscript above (full citation to be updated upon publication) and, where applicable, the original method references:

- Schneider, T., 2001: Analysis of incomplete climate data. *J. Climate*, 14, 853–871. https://doi.org/10.1175/1520-0442(2001)014%3C0853:AOICDE%3E2.0.CO;2
- Tierney, J. E., et al., 2015: Tropical sea surface temperatures for the past four centuries reconstructed from coral archives. *Paleoceanography*, 30, 226–252. https://doi.org/10.1002/2014PA002717

## License

Code in this repository is released under the license in [LICENSE](./LICENSE). Note that the RegEM dependency (linked above, not included here) is separately licensed under GPL-3.0.

## Contact

Yunfan Chen — yunfanch@umich.edu
