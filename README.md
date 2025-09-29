<h2>pathwayPCM</h2>

### Sample Codes

This repository contains an example input file in the `examples` directory so users can quickly try reconstructing macro-evolutionary pathway of multiple traits using pathwayPCM step-by-step:

**Step 1: Dataset Generation**

Input:

[`example.tree`](https://github.com/IwasakiLab/Evodictor/tree/main/example/example.tree): A phylogenetic tree in a Newick format.

[`OG_node_state.txt`](https://github.com/IwasakiLab/Evodictor/tree/main/example/OG_node_state.txt): The presence/absence profile of every ortholog group (OG) for every tip node (extant species) and every internal node (ancestors) of [`example.tree`](https://github.com/IwasakiLab/Evodictor/tree/main/example/example.tree). There is one row for every internal/tip node in this file. The first, second, and third columns of every row indicate the OG name, node name, and the presence/absence state, respectively. The presence/absence state is represented as `0` (absent), `1` (present), or `0.5` (uncertain; for ancestors). Rows for which states are `0` can be omitted in this file (in other words, states of nodes not defined in this file are treated as `0`).

[`example.tree`](https://github.com/IwasakiLab/Evodictor/tree/main/example/example.tree): Character coding table

Output:

[`branch_X_y.txt`](https://github.com/IwasakiLab/Evodictor/tree/main/example/output/branch_X_y.txt): The dataset for machine learning which can be an input file of `evodictor predict`. The first row is the header, and each of the following rows correspond to a branch in the [`example.tree`](https://github.com/IwasakiLab/Evodictor/tree/main/example/example.tree). The first, second, and third column of every row indicate the node name of a parental species of a branch in [`example.tree`](https://github.com/IwasakiLab/Evodictor/tree/main/example/example.tree), the number of present traits of every feature in the parental species (separated by `;`), and the occurrence of gene gain of predicted OG ([K00005](https://www.genome.jp/dbget-bin/www_bget?ko:K00005)) at the branch (`1`: the gene was gained at the branch;  `0`: the gene was not gained at the branch). 

**Step 2: Preparation of a command file for BayesTraits analysis**

Select top-20 important input features based on ANOVA F-value to predict gene gain of an OG ([K00005](https://www.genome.jp/dbget-bin/www_bget?ko:K00005)).

```shell
evodictor select -i branch_X_y.txt --skip_header --o1 feature_importance.txt --o2 selection_result.txt --o3 branch_X_y.selected.txt -k 20
```

