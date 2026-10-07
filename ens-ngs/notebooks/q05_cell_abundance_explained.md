### Why do we use the `q05` cell abundance?

Cell2location is a **Bayesian** model. For every spot and every cell state it does not return a single number, but a whole **posterior distribution** of plausible abundances (we draw 1000 samples from it with `export_posterior`). From this distribution the package stores several summaries, each with its own prefix:

| Prefix | Meaning |
|---|---|
| `means_` | mean of the posterior distribution |
| `stds_` | standard deviation |
| `q05_` | **5% quantile**: the abundance is above this value with ~95% probability |
| `q95_` | 95% quantile |

We work with `q05`, a **conservative** estimate. If the model is uncertain about a cell state in a spot, the distribution is wide and `q05` is pushed towards zero; if the model is confident, `q05` stays close to the mean. A high `q05` therefore means *"abundant and well supported by the data"*.

**This choice follows the original paper.** In the Methods of Liu et al. (*Spatiotemporal single-cell roadmap of human skin wound healing*, Cell Stem Cell 2025), section *"Deconvolution of ST-seq data and wound bulk RNA-seq data"*, the authors write:

> "The posterior distribution of cell abundance for each cell type in each spot was summarized as 5% quantile, representing high confidence, which was used for visualization and colocalization analysis."

The same section reports the two parameters that were changed from the defaults, which we also use:

- `N_cells_per_location = 20`: expected number of cells per spot, estimated by the authors from the average number of nuclei per spot in the H&E images.
- `detection_alpha = 20`: regularisation of the per-location normalisation, to account for large differences in RNA detection sensitivity between Visium spots.

#### Where the information comes from

- **From the paper (verified in its Methods):** the use of Cell2location, the 5% quantile as "high confidence" summary, `N_cells_per_location=20`, `detection_alpha=20`.
- **General knowledge about Cell2location (not stated in this paper):** the meaning of the posterior summaries (`means_`, `stds_`, `q05_`, `q95_`) and the interpretation of `q05` as a conservative estimate. See Kleshchevnikov et al., *Nature Biotechnology* 2022.

#### Differences between this notebook and the paper

This notebook follows the paper's Cell2location settings, but it is **not an exact reproduction**:

- The paper trained the reference regression model on **each cell type's top 100 marker genes**; here we select genes with Cell2location's `filter_genes` function instead.
- Our reference has **17 cell states** derived from our own clustering (`cell_state`), not the paper's annotation.
- Depending on the run, we may have used fewer training epochs and a subsampled reference for speed.

Expect similar *patterns*, not identical numbers.

*Question for you:* pick one cell state and plot both its `means_` and `q05_` abundance on the tissue. Where do the two maps differ most, and what does that tell you about the model's uncertainty?
