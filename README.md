# DOM-in-MicroEco

Dissolved organic matter (DOM), which is widely distributed in complex natural environments, plays a crucial role in global biogeochemical processes. Fourier transform ion cyclotron resonance mass spectrometry (FT-ICR-MS) is routinely applied to characterize the molecular composition of DOM. [Formularity](https://pnnl-comp-mass-spec.github.io/Formularity/) is one of many programs that have been developed for assigning molecular formulae in complex mass spectra.

Here, [FtmsAnalysis](https://github.com/KaiMa-endeavour/DOM-in-MicroEco/releases/tag/FtmsAnalysis) was developed to integrate molecular information from multiple samples. Please click [here](https://github.com/KaiMa-endeavour/MicroEcoTk/wiki) for more details.

## Dependency on MicroEcoTk

This project relies on the [**MicroEcoTk**](https://github.com/KaiMa-endeavour/MicroEcoTk) R package for common microbial ecology functions. Install it before running the scripts:

```r
install.packages("remotes")
remotes::install_github("KaiMa-endeavour/MicroEcoTk")
```

### Function Migration

The following functions were previously defined locally in this project's `Rscripts/` directory but have been moved to the MicroEcoTk package. They are no longer maintained here — please use the package versions instead.

| Former local function | MicroEcoTk equivalent | Description |
|---|---|---|
| `bNTI()` | `MicroEcoTk::bNTI()` | Beta nearest taxon index |
| `cohesion()` | `MicroEcoTk::cohesion()` | Community positive/negative cohesion |
| `commDispersal()` | `MicroEcoTk::commDispersal()` | Metacommunity dispersal index |
| `commNicheWidth()` | `MicroEcoTk::commNicheWidth()` | Levins' niche width per sample |
| `dominance()` | `MicroEcoTk::dominance()` | Dominance index |
| `fitsad_zipf()` | `MicroEcoTk::fitsad_zipf()` | Fit Zipf distribution to SAD |
| `min_error()` | `MicroEcoTk::min_error()` | Find scaling multiple with smallest beta-error |
| `normScale()` | `MicroEcoTk::normScale()` | Rescale vector to a target range |
| `nullModel()` | `MicroEcoTk::nullModel()` | Generate null community dissimilarity |
| `nullStoc()` | `MicroEcoTk::nullStoc()` | Community stochasticity via null models |
| `rank_abun()` | `MicroEcoTk::rank_abun()` | Rank abundance vector |
| `rarefy()` | `MicroEcoTk::rarefy()` | Rarefy community table to uniform depth |
| `rarefy_vt()` | `MicroEcoTk::rarefy_vt()` | Rarefy a single abundance vector |
| `rel_ab()` | `MicroEcoTk::rel_ab()` | Relative abundance (percentage) |
| `rarity()` | `MicroEcoTk::rarity()` | Skewness (rarity) of abundance vector |
| `Rsquare()` | `MicroEcoTk::Rsquare()` | Coefficient of determination (log10 scale) |
| `se()` | `MicroEcoTk::se()` | Standard error of the mean |
| `summStoc()` | `MicroEcoTk::summStoc()` | Summarise stochasticity indices |
| `taxa_partition()` | `MicroEcoTk::taxa_partition()` | Partition community into abundant/rare/occasional |

### Scripts that now use MicroEcoTk

These scripts previously contained local copies of the above functions. They now load `library(MicroEcoTk)` and call the package versions directly:

| Script | Package functions used |
|---|---|
| `SAD_models.R` | `Rsquare()`, `se()`, `rank_abun()`, `fitsad_zipf()` |
| `normScale.R` | `normScale()`, `min_error()` |
| `RandomSample_alphaDiversity.R` | `rarefy()`, `dominance()`, `rarity()`, `alpha_diversity()` |
| `keyMole.R` | `normScale()` |
| `DDR.R` | `taxa_partition()`, `rel_ab()` |
| `LDG.R` | `taxa_partition()`, `alpha_diversity()` |

### Quick usage examples

```r
library(MicroEcoTk)

# --- Rarefaction ---
rare_table <- rarefy(otu_table, depth = 10000)

# --- Alpha diversity ---
result <- alpha_diversity(comm, sample_names = rownames(comm),
                          methods = c("Richness", "Shannon", "Simpson"))

# --- Taxa partitioning ---
partition <- taxa_partition(comm, mode = "frequency")

# --- Null model stochasticity ---
stoc <- nullStoc(comm, null_model = "region", reps = 1000, nworker = 4)

# --- bNTI ---
result <- bNTI(comm, tree, reps = 1000, nworker = 4)

# --- Cohesion ---
result <- cohesion(comm, method = "spearman", nworker = 4)

# --- Niche width ---
nw <- commNicheWidth(comm, nworker = 4)

# --- Community dispersal ---
dispersal <- commDispersal(comm, d = geo_dist_matrix, nworker = 4)

# --- Dominance & rarity ---
d <- dominance(abundance_vector)
r <- rarity(abundance_vector)

# --- Scaling for DOM intensity ---
scaled <- normScale(intensity_vector, floor = 100, multiple = 9.8e-5)
```

For full documentation of each function, see the MicroEcoTk package help pages (e.g., `?MicroEcoTk::bNTI`) or visit the [package repository](https://github.com/KaiMa-endeavour/MicroEcoTk).

## References

Ma, K., Li, Y., Song, W., Zhou, J., Liu, X., Wang, M., Gong, X., Wang, L., Tu, Q., 2024. Disentangling drivers of mudflat intertidal DOM chemodiversity using ecological models. Nat. Commun. 15, 6620. https://doi.org/10.1038/s41467-024-50841-9

Ma, K., Li, Y., Liu, X., Song, W., Zhou, J., Gong, X., Wang, M., Li, C., Liu, J., Tu, Q., 2023. Bacteria rather than fungi mediate the chemodiversity of dissolved organic matter in a mudflat intertidal zone. Science of The Total Environment 893, 164835. https://doi.org/10.1016/j.scitotenv.2023.164835
