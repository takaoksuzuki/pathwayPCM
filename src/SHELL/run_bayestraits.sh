#!/bin/bash

# Run BayesTraits for one phylum
./BayesTraits/BayesTraitsV4 \
  data/treefile.tree \
  data/traitfile.tsv \
  < results/inputcommands_Proteobacteria.txt \
  > results/output_Proteobacteria.txt
