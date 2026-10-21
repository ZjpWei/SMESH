# SMESH analysis

## Introduction

This project includes the analyses of `SMESH` paper. Please follow the `README.md` file in each folder to reproduce the simulation and three real data applications. 
The analyses use `miMeta` version `2.1.0` and `PALM` version `0.1.1`. Previous versions are available on the [miMeta](https://github.com/ZjpWei/miMeta/releases) and [PALM](https://github.com/ZjpWei/PALM_package/releases) release pages.

SMESH is an end-to-end pipeline that combines reliable within-study estimation of AA association effects and their uncertainty with a 
new structure-learning meta-analysis model. Unlike methods that impose a common effect, require prespecified context groups, or separate 
clustering from meta-analytic signature selection, SMESH jointly learns association-defined context clusters and sparse cluster-level effects.

The repository holds three real-data applications and two simulation studies:

| Application | Contexts | Taxa |
|---|---|---|
| **Pan-disease** | 12 diseases, 32 datasets | 244 species |
| **Pan-tumor** | 13 tumour types, 28 datasets | 157 genera |
| **Colorectal neoplasia** | 26 study × stage strata | 260 species |

---

## Structure

```
SMESH/
    ├── Analysis/           the pipeline        — see README file for more details
    ├── Rscript/            figures and tables  — see README file for more details
    ├── utility/            the utility files
    ├── Data/               raw data
    ├── Analysis/           the pipeline        — run these, in order
    ├── Rscript/            figures and tables  — run these afterwards
    ├── utility/            the estimator and shared helpers
    ├── Data/               raw inputs, exactly as downloaded
    ├── GMrepo_analysis/    pan-disease  — processed data, fits, panels
    ├── HGMT_analysis/      pan-tumor    — processed data, fits, panels
    ├── CRC_analysis/       colorectal neoplasia   — processed data, fits, panels
    ├── Simulation/         simulation template, replicates and panels
    └── Figure/             the assembled figures and supplementary tables
```

## References

* Wei Z, Chen G, Tang ZZ. Structure-learning meta-analysis resolves cross-context heterogeneity in microbiome associations.

## Contact

If you have questions about the code or the analyses presented in this project, please feel free to contact me [Weizhou](mailto:zwei74@wisc.edu?subject=[GitHub]%20Melody%20paper%20analysis).

