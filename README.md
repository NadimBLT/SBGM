# SBGM
Stratified Binary Graphical Models

## Packages required 


```
MatrixModels, glmnet, BMN and FusedLasso.
```
Note that the FusedLasso package is not maintained on the CRAN anymore. You can download it from source file which is attached in this file, but it work only on linux.



## Usage
### SepLogit_Fused, SepLogit_DataShared, SepLogit_Red and SepLogit_Indep
#### Arguments
* **data** : a matrix whose last column corresponds to a categorical variable that defines the strata.
* **adap** : adapt = FALSE corresponds to the standard version of the methods and adapt = TRUE corresponds to the adaptive version of the method.
* **standardize** : Default is standardize=FALSE. If standardize=FALSE, we standardize the variables.

#### Value
* **MIN** :  graphical models returned by the SepLogit_MIN and BIC critera.
* **MAX** :  graphical models returned by the SepLogit_MAX and BIC critera.

### Function_Guo
#### Arguments
* **data** : a matrix whose last column corresponds to a categorical variable that defines the strata.

#### Value
* **BIC** :  graphical models returned by BIC criterion.

