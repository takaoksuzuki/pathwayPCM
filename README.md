<h1>pathwayPCM</h1>

## Overview of pathwayPCM

This method models trait evolution as a continuous-time Markov process over a state space of trait combinations, where nodes correspond to trait combination codes and directed edges denote evolutionarily permissible transitions between them.

<img src=image/fig1.jpg >

(A) Model constraints. Nodes represent all possible combinations of three traits, corresponding to trait combination codes (here illustrated with binary traits), and edges represent evolutionarily possible transitions between them. The initial model assumes a fully connected topology with constrained rate parameters—six transition rates (three forward, three backward) indicated by distinct colors. Although illustrated here for three binary traits, the pathwayPCM framework can also handle multistate (non-binary) traits. 
(B) Bayesian sampling by reversible-jump MCMC (rjMCMC). The rjMCMC algorithm simultaneously updates transition rates (edge thickness) and network topology (edge addition/deletion). After convergence, rate parameters are averaged, and the topology is determined based on z-scores of edge inclusion frequencies, yielding an optimized macroevolutionary transition network (the macroevolutionary pathway). 
(C) Downstream analysis. Because the inferred macroevolutionary pathway has a network structure, standard network analyses can be applied, such as module structure detection, to identify groups of traits that tend to evolve together. Example modules are shown in pink and yellow.

<img width="13361" height="53" alt="image" src="https://github.com/user-attachments/assets/0699977f-ff7f-48cb-b99c-28e0e2dfd5c4" />

## System requirements

- R versions ≥ 4.4.
- Required R packages:
  - readr (≥ 2.1.5)
  - dplyr (≥ 1.1.4)
  - tidyr (≥ 1.3.1)
  - stringr (≥ 1.5.1)
  - igraph (ver. 2.0.3)
- Tested on Linux (Ubuntu 22.04) and Windows 11
- BayesTraits version 4.0 or later (command-line executable accessible in $PATH)

## pathwayPCM Workflow

pathwayPCM consists of three computational stages:

> [!NOTE]
> pathwayPCM is currently implemented as a modular pipeline, rather than a single integrated platform.
> This is because the method relies on BayesTraits v4, which must be executed as an external command outside R.
> Each step is therefore executed separately (via R scripts and shell scripts). A unified interface is planned for future releases.

#### 1. Data preprocessing
(corresponds to Sample Codes Steps 1–2 using R and shell scripts)

#### 2. Phylogenetic comparative analysis using reversible-jump MCMC
(corresponds to Sample Codes Steps 3–4 using BayesTraits and shell scripts)

#### 3. Downstream analysis
(corresponds to Sample Codes Steps 5–6 using R, including igraph’s Infomap implementation)

## Installation Guide

**Typical install time on a "normal" desktop computer:** ~5 minutes.

#### 1. Prerequisites
pathwayPCM is currently implemented as a modular pipeline, not a single integrated R package.

To execute the workflow, users must install:
- R (≥ 4.4)
- Required R packages
- BayesTraits v4 (external executable)
- (Optional) Shell environment (bash)
Once installed, users can run the three computational steps (data preprocessing → PCMs → downstream analysis).

#### 2. Installation instructions
##### Step 1 — Clone the repository
```SHELL
git clone https://github.com/takaoksuzuki/pathwayPCM.git
cd pathwayPCM
```
This repository contains the R scripts, shell scripts, and example input files necessary to reproduce the workflow.

##### Step 2 — Install R and required R packages
Install R (≥ 4.4) from:
https://cran.r-project.org/
Then install required packages:
```R
install.packages(c(
  "readr", 
  "dplyr", 
  "tidyr", 
  "stringr", 
  "igraph"
))
```
These packages are used in the preprocessing and downstream analysis steps (Sample Codes Steps 1, 4–6).

##### Step 3 — Prepare shell environment (bash)
Some steps in the pathwayPCM workflow (especially running BayesTraits) require a working bash environment.

To check that bash is available:
```SHELL
bash --version
```
Linux and macOS provide a bash shell by default. Windows users should run the pipeline through WSL2 for full compatibility.

##### Step 4 — Install BayesTraits v4
BayesTraits v4 must be downloaded separately because it is not an R package.

Download (free academic license):
https://www.evolution.reading.ac.uk/BayesTraits.html

Unzip the downloaded archive and place the BayesTraitsV4 executable in:
```SHELL
mv /path/to/downloaded/BayesTraitsV4 ./BayesTraits/
chmod +x ./BayesTraits/BayesTraitsV4
```
Test the installation:
```SHELL
./BayesTraits/BayesTraitsV4
```
BayesTraits is used in Sample Codes Step 3 (“Run BayesTraits”) of the pipeline.

## Sample Codes

**Expected run time for demo on a "normal" desktop computer:** ~10–20 seconds.
Full BayesTraits analysis requires an HPC system (~7 days, 1 CPU).

This repository contains an example input file in the `data` directory so users can quickly try reconstructing macro-evolutionary pathway of multiple traits using pathwayPCM step-by-step:

**Step 1: Dataset Preparation**

Input:

[`Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/data/Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree):<BR>
A phylogenetic tree in a Nexus format.

[`ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/data/ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv):<BR>
Input table of accession–trait combinations. The first column lists tip labels, and the second column contains alphanumeric character codes (e.g., A, B, C, 1, 2, …) that uniquely identify combinations of traits (cell shape, motility, and sporulation).

[`ch_coding_GTDB207_cellshape_motility_sporulation_ALLphylum.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/data/ch_coding_GTDB207_cellshape_motility_sporulation_ALLphylum.tsv):<BR>
Character coding table. This file defines the correspondence between each alphanumeric code and the actual combination of trait states (e.g., A = rod–motile–sporulating).
Used as a unified reference across phyla.

**Step 2: Prepare BayesTraits command files (w/ rate restrictions)**

`makerestriction.R` generates intermediate restriction files (mycommand_*.txt) for the following step. These files are saved in the [`results/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/results) directory.

```R
./src/R/makerestriction.R
```
`makeinputcommands.R` creates the BayesTraits input command files for Step 3 by integrating the restriction files. Output files are saved as: inputcommands_*.txt in the [`results/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/results) directory.
```R
./src/R/makeinputcommands.R
```

**Step 3: Run BayesTraits**

Preparing shell files.

```SHELL
./src/SHELL/run_bayestraits.sh
```

> [!CAUTION]
> Running BayesTraits on the full dataset generates >13 GB of output
> and requires approximately 7 days of computation on an HPC cluster (1 CPU).
> We strongly recommend running this step on a high-performance computing system.

**Step 4: Parse the results of BayesTraits**

In batch processing mode, each script (`extractdata.R`, `makeqrate.R`, `makezscore.R`) automatically loads the BayesTraits output files located in the [`results/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/results) directory, which are generated when users run Step 3 (run_bayestraits.sh).

Parse the BayesTraits output files to extract data rows.
```R
./src/R/extractdata.R
```
Calcurate qrates from the extracted data.

```R
./src/R/makeqrate.R
```
Calcurate z-scores from the extracted data.

```R
./src/R/makezscore.R
```

> [!Note]
> If you prefer to skip the long BayesTraits computation, pre-computed example results
> are provided in the `output/` directory.  
> Running Step 4 and Step 5 using these files allows you to reproduce the downstream
> analyses without executing Step 3 (`run_bayestraits.sh`).


Output:

If you use the pre-computed example results in the `output/` directory (instead of running Step 3),  
the following files are available:

- [`qrate_ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/output/qrate_ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv):<BR>

  Estimated transition rates (q-values) between all multi-trait states.

- [`zscore_ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/output/zscore_ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv):<BR>
  Percentage of MCMC samples in which each transition rate was constrained to zero  
  (i.e., number of zero-constrained iterations / post-burn-in iterations × 100).

**Step 5: Draw the macro-evolutionary pathway**

This script reads the q-rate and z-score tables from [`results/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/results) (or [`output/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/output) if Step 3 is skipped) and draws the macro-evolutionary pathway network. Output figures are saved in `results/`.

```R
./src/R/plot_zscore+qrates.R
```

Output:

<img src=image/fig2.jpg >

**Step 6: Detecte the flow modules by Infomap methods**

Identify flow modules (clusters of trait transitions) within the macro-evolutionary pathway network using the Infomap algorithm.

```R
./src/R/calc_module_infomap.R
```

Output

Example output files generated by `calc_module_infomap.R` are provided in the [`output/`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/output) directory

[`modulePathway_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/output/modulePathway_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv):<BR>
Table listing the Infomap-detected module assignments for each trait-state code.  
The first column contains alphanumeric state codes (e.g., A, B, C, 1, 2, …), and the second column indicates the module ID to which each state belongs.

## Licence

The pathwayPCM code is released under the MIT License.
