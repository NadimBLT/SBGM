# SBGM
Stratified Binary Graphical Models

## Packages required 


```
MatrixModels, glmnet, BMN and FusedLasso.
```
Note that the FusedLasso package is not maintained on the CRAN anymore. You can download it from source file which is attached in this file, but it work only on linux.



## Usage
### SepLogit_Fused, SepLogit_DataShared, SepLogit_Ref and SepLogit_Indep
#### Arguments
* **data**        : input matrix, of dimension n x (p + 1). n is the number of observations and (p + 1) is the number of variables, where                     the first p columns correspond to p binary variables and the last column corresponds to a categorical variable                           defining the strata.
* **adap**        : If True, the L1-norm penalty terms are replaced by weighted terms. The weights are based on estimates obtained by                         applying OLS method independently for each stratum. Default is adap=FALSE.
* **standardize** : If True, the original variables are standardized. Default is standardize=FALSE.

#### Value
* **MIN** :  graphical models returned by the SepLogit_MIN and BIC critera.
* **MAX** :  graphical models returned by the SepLogit_MAX and BIC critera.

### Function_Guo
#### Arguments
* **data** : input matrix, of dimension n x (p + 1). n is the number of observations and (p + 1) n is the number of variables, where the first p columns correspond to the p binary variables and the last column corresponds to a categorical variable defining the strata.

#### Value
* **BIC** :  graphical models returned by BIC criterion.

