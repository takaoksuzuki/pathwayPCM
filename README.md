<h1>pathwayPCM</h1>

## Overview of pathwayPCM

This method models trait evolution as a continuous-time Markov process over a state space of trait combinations, where nodes correspond to trait combination codes and directed edges denote evolutionarily permissible transitions between them.

<img src=image/fig1.jpg >

## System requirements

This method models trait evolution as a continuous-time Markov process over a state space of trait combinations, where nodes correspond to trait combination codes and directed edges denote evolutionarily permissible transitions between them.

## Sample Codes

This repository contains an example input file in the `Example` directory so users can quickly try reconstructing macro-evolutionary pathway of multiple traits using pathwayPCM step-by-step:

**Step 1: Dataset Preparation**

Input:

[`Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/Example/Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree): A phylogenetic tree in a Nexus format.

[`ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/Example/ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv): Character coding table

**Step 2: Preparation of a command file for BayesTraits analysis**

Make command lines for constraining transition rates.

```R
/src/R/makerestriction.R
```

**Step 3: Run BayesTraits**

Preparing shell files.

```SHELL
/src/SHELL/run_bayestraits.sh
```

**Step 4: Parse the results of BayesTraits**

Parsing the result files to extract data rows. Because of batch processing, users do not need to specify input files.

```R
/src/R/extractdata.R
```
Calcurate qrates from the extracted data.

```R
/src/R/makeqrate.R
```
Calcurate z-scores from the extracted data.

```R
/src/R/makezscore.R
```

**Step 5: Draw the macro-evolutionary pathway**

Make edge list of macro-evolutioary pathway and draw it.

```R
/src/R/plot_zscore+qrates.R
```


