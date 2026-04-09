# myopia-ml

# Myopia Chronobiological Signature ML Framework

This repository contains the data and R scripts used to discover and validate a chronobiological gene signature for myopia, as described in preprint doi: https://doi.org/10.64898/2026.04.02.716020.

## Repository Structure
* `myopia_model_stages1_to_3.R`: Script for Stage 1 (Retina Onset Discovery), Stage 2 (Choroid Cross-Tissue Validation), and Stage 3 (Progression Cross-Stage Validation).
* `external_validation_stage4.R`: Script for Stage 4 (Independent External Validation).
* `my_53_gene_signature.txt`: The final consensus genes selected via Boruta/LASSO.
* `*.csv`: Raw count matrices and phenotype metadata for the various datasets.

## Requirements
To run these scripts, you will need R and the following packages:
`edgeR`, `tidyverse`, `caret`, `randomForest`, `pROC`, `glmnet`, `Boruta`, `e1071`.

## Usage
1. Download or clone this repository.
2. Open `myopia_model_stages1_to_3.R` and run to generate the primary models and inner validations. 
3. Run `external_validation_stage4.R` using the generated signature to validate against the external GSE203604 dataset.
