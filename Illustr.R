path_repository="~/SBGM-master/"


source(paste0(path_repository,"Functions-For-SepLogit.R"))
source(paste0(path_repository,"SepLogit-Ref-DataShared.R"))
source(paste0(path_repository,"SepLogit-Indep.R"))
source(paste0(path_repository,"Functions_For_Fused.R"))
source(paste0(path_repository,"SepLogit-Fused.R"))
source(paste0(path_repository,"Fonctions_For_Guo.R"))
source(paste0(path_repository,"Function_Guo.R"))

#install.packages(paste0(path_repository,"FusedLasso_1.0.6.tar.gz"),repos=NULL,type="source")
#install.packages(paste0(path_repository,"BMN_1.02.tar.gz"),repos=NULL,type="source")



library(MatrixModels)
library(glmnet)
#library(FusedLasso)
#library(BMN)


load(paste0(path_repository,"/M-3NN-R-0.25.Rdata"))
Data=fichier[[2]]
head(Data)



#RES_GUO          =Function_Guo(data=Data)
#RES_FUSED        =SepLogit_Fused(data = Data,adap = FALSE,standardize = FALSE)
RES_INDEP        =SepLogit_Indep(data = Data,adap = FALSE,standardize = FALSE)
RES_REF          =SepLogit_Ref(data = Data,adap = FALSE,standardize = FALSE)
RES_DataShared   =SepLogit_DataShared(data = Data,adap = FALSE,standardize = FALSE)
