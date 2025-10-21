<h2>pathwayPCM</h2>

### Overview of pathwayPCM

This method models trait evolution as a continuous-time Markov process over a state space of trait combinations, where nodes correspond to trait combination codes and directed edges denote evolutionarily permissible transitions between them.

<img src=image/fig1.jpg >

### Sample Codes

This repository contains an example input file in the `Example` directory so users can quickly try reconstructing macro-evolutionary pathway of multiple traits using pathwayPCM step-by-step:

**Step 1: Dataset Preparation**

Input:

[`Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/Example/Proteobacteria_cellshape_motility_sporulation_bac120_r207.tree): A phylogenetic tree in a Nexus format.

[`ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv`](https://github.com/takaoksuzuki/pathwayPCM/blob/main/Example/ch_GTDB207_cellshape_motility_sporulation_Proteobacteria.tsv): Character coding table

**Step 2: Preparation of a command file for BayesTraits analysis**

Make command lines for constraining transition rates.

```R
makerestriction.R
```

**Step 3: Run BayesTraits**

Preparing shell files.

```SHELL
bac_16.sh
```

**Step 4: Parse the results of BayesTraits**

Parsing the result files to extract data rows. Because of batch processing, users do not need to specify input files.

```R
extractdata.R
```
Calcurate qrates from the extracted data.

```R
makeqrate.R
```
Calcurate z-scores from the extracted data.

```R
makezscore.R
```

**Step 5: Draw the macro-evolutionary pathway**

Make edge list of macro-evolutioary pathway and draw it.

```R
evodictor select -i branch_X_y.txt --skip_header --o1 feature_importance.txt --o2 selection_result.txt --o3 branch_X_y.selected.txt -k 20
```


