# Kernel PCA for Multi-Omics Data Integration

This repository provides code for performing kernel PCA and association testing on multi-omics data (e.g., genotypes and methylation) using kernel machine regression.

## Prerequisites

You will need the following R packages:

    - GENESIS: For fitting null models and handling related individuals.
    - SKAT: For kernel-based association testing.
    - CompQuadForm: For Davies' method and saddlepoint approximation.

Install these packages using the following R command:

``` r
install.packages(c("GENESIS", "SKAT", "CompQuadForm"))
```

## File Structure

- `kernel_analysis.R`: The main script to run kernel PCA analysis and association testing.
- `kernel_helpers.R`: Contains helper functions to calculate test statistics, p-values, and kernel matrix operations.
- `data/`: Directory to store input genotype, methylation, and phenotype data files (user provided).

## How to Use

1. **Clone the repository:**

``` bash
git clone https://github.dev/amarise/omics_association_testing/
```

2. **Prepare your data:**

    - Input files should include:
        - Genotype matrix (`G`)
        - Methylation matrix (`M`)
        - Outcome variable (`y`)
        - Covariates (`X`)
        - Genetic relatedness matrix (`grm`)
    - Ensure that the row order in each data file is conserved, so they correspond to the same subjects across all files.

3. **Run the analysis:**

    - Source the kernel_analysis.R script in R:

``` r
source("kernel_analysis.R")
```

    - Modify the script to point to your data files, then run the kernel PCA analysis:

``` r
result <- kernelPCA(w = seq(0, 1, by = 0.1), y = y_data, X = covariates, G = genotype_data, M = methylation_data, grm = relatedness_matrix)
```
4. **Review Results:**
    
    - The results of the analysis will include p-values, test statistics, and the optimal weight (`w`) for combining kernels.

## Example Usage

Here’s an example of how to run the `kernelPCA` function:

``` r
# Example of running the analysis
result <- kernelPCA(w = seq(0, 1, by = 0.1), y = y_data, X = covariates, G = genotype_data, M = methylation_data, grm = relatedness_matrix)

# View the results
print(result$grid_results)  # View the grid of test statistics and p-values
print(result$p_val)         # Overall p-value for the test
print(result$w)             # Optimal weight for kernel combination
```

## Custom Functions

    - `kernelPCA()`: This function performs kernel PCA and association testing for two omics data types (genotypes and methylation).
    - `SKAT_2Kernel_Optimal_Get_Q()`: Computes the test statistics for a range of kernel weights.
    - `SKAT_2Kernel_Ortho_Optimal_Get_Pvalue()`: Calculates p-values for kernel association tests.

## Acknowledgments

This code uses methods from the [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) and [SKAT](https://cran.r-project.org/web/packages/SKAT/index.html) R packages. Please cite relevant papers if you use this code in your research. ​
​