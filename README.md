# SBGM
Stratified Binary Graphical Models

Functions SepLogit_DataShared (available in file SepLogit-Ref-DataShared.R), SepLogit_Fused (available in file SepLogit_Fused.R),  SepLogit_Ref (available in file SepLogit-Ref-DataShared.R), SepLogit_Indep (available in file SepLogit_Indep) and Guo (available in file Function_Guo.R) are used to estimate several binary graphical models defined on several strata.

SepLogit_Indep is based on the SepLogit approach and estimates the models in an independent way by applying the seplogit approach on each stratum separately.

Guo is based on the pseudo like-lihood method and estimates the models in a semi-joint way using a multiplicative decomposition of parameters.

SepLogit_Ref is based on the seplogit approach and estimates the models in a semi-joint way by choosing a reference stratum a priori and by using an additive decomposition of parameters.

SepLogit_Fused is based on the seplogit approach and estimates the models in a joint way using the fused penalty.

SepLogit_Datashared is based on the seplogit approach and estimates the models in a joint way using the DataShared method.

See Illustr.R for a simple example illustrating the use of these functions.

See https://arxiv.org/pdf/1709.10298.pdf for more details.
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
* **MIN**         : graphical models returned by the SepLogit_MIN and BIC critera.
* **MAX**         : graphical models returned by the SepLogit_MAX and BIC critera.

### Function_Guo
#### Arguments
* **data**        : input matrix, of dimension n x (p + 1). n is the number of observations and (p + 1) n is the number of variables,                         where the first p columns correspond to the p binary variables and the last column corresponds to a categorical                           variable defining the strata.

#### Value
* **BIC**         : graphical models returned by BIC criterion.

