# Association Testing for Two Multi-Omics Data Types

This repository provides R code for performing association testing on two multi-omics data types (e.g., genotypes and methylation) using two methods: kernel PCA and perturbation testing. These methods allow for the assessment of associations between multi-omics data and a trait of interest.

## Prerequisites

The following R packages are required:

- `GENESIS`: For fitting null models and handling related individuals.
- `SKAT`: For kernel-based association testing.
- `CompQuadForm`: For Davies' method and saddlepoint approximation.

Install the required packages using the command below:

``` r
install.packages(c("GENESIS", "SKAT", "CompQuadForm"))
```

## File Structure

- `kernelPCA_method.R`:  Script for running the kernel PCA association testing method.
- `perturbation_method.R`: Script for running the perturbation method for association testing.
- `kernel_helpers.R`: Contains helper functions to calculate test statistics, p-values, and kernel matrix operations.
- `data/`: Directory to store input omics data type 1 (e.g. genotype), omics data type 2 (e.g. methylation), and phenotype data files (user provided).

## How to Use

1. **Clone the repository:**

``` bash
git clone https://github.dev/amarise/omics_association_testing/
```

2. **Prepare your data:**

Ensure that the following data files are prepared and correspond to the same subjects:
- Omics data type 1 matrix (`omics1`)
- Omics data type 2 matrix (`omics2`)
- Outcome variable (`y`)
- Covariates (`X`)
- Genetic relatedness matrix (`grm`)
Make sure that the row order in each data file is consistent across all files, so that they correspond to the same subjects.

## Kernel PCA Association Testing

### Running Kernel PCA Analysis:

1. Source the `kernelPCA_method.R` script in R:

    ``` r
    source("kernelPCA_method.R")
    ```

2. Modify the script to point to your data files, then run the kernel PCA analysis:

    ``` r
    result <- kernelPCA(w = seq(0, 1, by = 0.25), y = y_data, X = covariates, omics1 = omics1_data, omics2 = omics2_data, grm = relatedness_matrix)
    ```

3. Review Results:
    
    - The results of the analysis will include p-values, test statistics, and the optimal weight (`w`) for combining kernels.

### Example Usage:

``` r
result <- kernelPCA(w = seq(0, 1, by = 0.25), y = y_data, X = covariates, omics1 = genotype_data, omics2 = methylation_data, grm = relatedness_matrix)

# View the results
print(result$grid_results)  # View the grid of test statistics and p-values
print(result$p_val)         # Overall p-value for the test
print(result$w)             # Optimal weight for kernel combination
```

## Perturbation Method for Association Testing

### Running Perturbation Method:

1. Include the kernelPCA_method.R file in your project.
2. Source the script and run the perturbation method using your data:

```r
source("perturbation_method.R")
result <- perturbation(w = seq(0, 1, by = 0.25), y = y_data, X = covariates, omics1 = omics1_data, omics2 = omics2_data, grm = relatedness_matrix)
```

3. Review Results:
    - This method will provide test statistics, p-values, and an optimal weight (`w`) for combining kernels based on perturbations.

### Example Usage:

```r
result <- perturbation(w = seq(0, 1, by = 0.1), y = y_data, X = covariates, omics1 = genotype_data, omics2 = methylation_data, grm = relatedness_matrix)

# View the results
print(result$grid_results)  # View the test statistics and p-values for different weights
print(result$p_val)         # Overall p-value for the test
print(result$w)             # Optimal weight for kernel combination
```

## Custom Functions

- `kernelPCA()`: Performs kernel PCA and association testing for two omics data types.
- `perturbation()`: Implements the perturbation-based association testing method.
- `SKAT_2Kernel_Optimal_Get_Q()`: Computes test statistics for a range of kernel weights.
- `SKAT_2Kernel_Ortho_Optimal_Get_Pvalue()`: Calculates p-values for kernel association tests.

## Acknowledgments

This code uses methods from the [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) and [SKAT](https://cran.r-project.org/web/packages/SKAT/index.html) R packages. Please cite relevant papers if you use this code in your research. ​
​