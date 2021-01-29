screen.DataSharedLasso.GLMs.ITCPT = function(X, y, Z, p=ncol(X), N_Str=max(Z), family = "gaussian", Balanced= F, wObs= rep(1, length(y)), ratio = 0.001, Nlam = 50, wLambda1=NULL, wLambda2 = NULL)
{
  # X:        design matrix (n,p)
  # y:        vector of labels (n,1)
  # Z:        vector containing the stratum to which each observation belongs (n,1); must range from 1 to N_Str 
  # p:        number of predictors in X
  # N_Str:    number of strata (i.e., max(Z))
  # family:   "gaussian", "binomial", etc.
  # Balanced:	T or F: are the strata balanced or not.
  # wObs:		Observations weights
  # ratio:    lambda_min/lambda_max
  # Nlam: 	Lenght of the lambda sequences
  # wLambda1: adaptive weights for the common effects: vector of size p 
  # wLambda2: adaptive weights for the specific effects..Size Kp: 1,...p for the specific coef in the first stratum, p+1,...,2p : in the 2nd stratum
  
  
  library(glmnet)
  
  X         = X[order(Z),]
  y         = y[order(Z)]
  ncs       = as.numeric(table(Z)) 
  
  calX.Hier = Hierise_ITCPT(X, ncs)
  IND = which(calX.Hier!=0,arr.ind=T)
  calX.Hier = sparseMatrix(i = IND[,1], j = IND[,2] , x= calX.Hier[IND])
  
  tauks     = rep(1,N_Str)
  if (Balanced == F) {tauks = sqrt(N_Str)*sqrt(ncs/sum(ncs))} # so that tauks=1 in the balanced case
  
  if (is.null(wLambda1)) {wL1 = rep(1, p)}else {wL1 = wLambda1}
  if (is.null(wLambda2)) {wL2 = rep(tauks, each=p+1)}else {wL2 = wLambda2}
  
  wObs      = wObs*length(y)/sum(wObs)  # rescaling the weights so that they sum up to n
  
  #########################
  MatDiagWeightsInv = sparseMatrix(i = 1:(p*(N_Str+1)+N_Str), j = 1:(p*(N_Str+1)+N_Str) ,x= 1/c(wL1,wL2) )
  ratio.max       = 1.01*sum(tauks)
  
  ratio_lam1_lam2 = seq(ratio.max,ratio.max*ratio, length.out=Nlam) # Les résultats du papier ont été obtenu avec : make_lambda_seq(N_Str, ratio, Nlam), mais pas de raison de se mettre sur une echelle log; si ?
  #ratio_lam1_lam2=1
  ALL_COEFS = BIC = NULL
  for (rk in ratio_lam1_lam2)
  {
    rk_mat = sparseMatrix(i = 1:(p*(N_Str+1)+N_Str), j = 1:(p*(N_Str+1)+N_Str) , x=c(rep(1/rk,p), rep(1,N_Str*(p+1))) )
    calX.Hier.rat = calX.Hier%*%MatDiagWeightsInv%*%rk_mat
    
    mod0         = glmnet(calX.Hier.rat, y, weights= wObs, standardize=F, family=family, nlambda=Nlam, lambda.min.ratio=ratio)
    
    dev0         = mod0$nulldev*(1-mod0$dev.ratio)
    Devf         = dev0 - mod0$nulldev
    BIC          = c(BIC, Devf + mod0$df*log(length(y)))
    
    bbeta 		 = rbind(mod0$a0, as.matrix(mod0$beta))
    bbeta		 = Matrix(bbeta, sparse=T)
    
    MatDiagWeightsInv_ITCPT = Matrix(cbind(c(1, rep(0, nrow(MatDiagWeightsInv))), rbind(rep(0, ncol(MatDiagWeightsInv)), as.matrix(MatDiagWeightsInv))), sparse=T)
    
    ALL_COEFS = cbind(ALL_COEFS, apply(bbeta,2,function(v){mise_forme_coef_hiers(Matrix(diag(c(1, rep(1/rk,p), rep(1,N_Str*(p+1)))), sparse=T)%*%MatDiagWeightsInv_ITCPT%*%v,N_Str)}))
    
    
  }
  
  return(list(ALL_COEFS=ALL_COEFS, BIC=BIC))	
}






screen.RefLasso.GLMs.ITCPT = function(X, y, Z, p=ncol(X), N_Str=max(Z), family = "gaussian", Balanced= F, wObs= rep(1, length(y)), ratio = 0.001, Nlam = 50, wLambda1=NULL, wLambda2 = NULL)
{
  
  # X:        design matrix (n,p)
  # y:        vector of labels (n,1)
  # Z:        vector containing the stratum to which each observation belongs (n,1); must range from 1 to N_Str 
  # p:        number of predictors in X
  # N_Str:    number of strata (i.e., max(Z))
  # family:   "gaussian", "binomial", etc.
  # Balanced:	T or F: are the strata balanced or not.
  # wObs:		Observations weights
  # ratio:    lambda_min/lambda_max
  # Nlam: 	Lenght of the lambda sequences
  # wLambda1: adaptive weights for the common effects: vector of size p 
  # wLambda2: adaptive weights for the specific effects..Size Kp: 1,...p for the specific coef in the first stratum, p+1,...,2p : in the 2nd stratum
  
  
  library(glmnet)
  
  X         = X[order(Z),]
  y         = y[order(Z)]
  ncs       = as.numeric(table(Z)) 
  
  calX.Hier = Hierise_ITCPT(X, ncs)
  calX.Hier = calX.Hier[,-((p+1):(2*p+1))]
  calX.Hier=as.matrix(calX.Hier)
  IND = which(calX.Hier!=0,arr.ind=TRUE)
  calX.Hier = sparseMatrix(i = IND[,1], j = IND[,2] , x= calX.Hier[IND])
  N_Str0=N_Str
  N_Str=N_Str-1
  ncs0=ncs
  ncs=ncs[-1]
  tauks     = rep(1,N_Str)
  if (Balanced == F) {tauks = sqrt(N_Str0)*sqrt(ncs/sum(ncs0))} # so that tauks=1 in the balanced case
  
  if (is.null(wLambda1)) {wL1 = rep(1, p)}else {wL1 = wLambda1}
  if (is.null(wLambda2)) {wL2 = rep(tauks, each=p+1)}else {wL2 = wLambda2}
  
  wObs      = wObs*length(y)/sum(wObs)  # rescaling the weights so that they sum up to n
  
  #########################
  MatDiagWeightsInv = sparseMatrix(i = 1:(p*(N_Str+1)+N_Str), j = 1:(p*(N_Str+1)+N_Str) ,x= 1/c(wL1,wL2) )
  ratio.max       = 1.01*sum(tauks)
  
  ratio_lam1_lam2 = seq(ratio.max,ratio.max*ratio, length.out=Nlam) # Les résultats du papier ont été obtenu avec : make_lambda_seq(N_Str, ratio, Nlam), mais pas de raison de se mettre sur une echelle log; si ?
  #ratio_lam1_lam2=1
  ALL_COEFS = BIC = NULL
  for (rk in ratio_lam1_lam2)
  {
    rk_mat = sparseMatrix(i = 1:(p*(N_Str+1)+N_Str), j = 1:(p*(N_Str+1)+N_Str) , x=c(rep(1/rk,p), rep(1,N_Str*(p+1))) )
    calX.Hier.rat = calX.Hier%*%MatDiagWeightsInv%*%rk_mat
    
    mod0         = glmnet(calX.Hier.rat, y, weights= wObs, standardize=F, family=family, nlambda=Nlam, lambda.min.ratio=ratio)
    
    dev0         = mod0$nulldev*(1-mod0$dev.ratio)
    Devf         = dev0 - mod0$nulldev
    BIC          = c(BIC, Devf + mod0$df*log(length(y)))
    
    bbeta 		 = rbind(mod0$a0, as.matrix(mod0$beta))
    bbeta		 = Matrix(bbeta, sparse=TRUE)
    
    MatDiagWeightsInv_ITCPT = Matrix(cbind(c(1, rep(0, nrow(MatDiagWeightsInv))), rbind(rep(0, ncol(MatDiagWeightsInv)), as.matrix(MatDiagWeightsInv))), sparse=TRUE)
    
    ALL_COEFS = cbind(ALL_COEFS, apply(bbeta,2,function(v){mise_forme_coef_hiers_R(Matrix(diag(c(1, rep(1/rk,p), rep(1,N_Str*(p+1)))), sparse=T)%*%MatDiagWeightsInv_ITCPT%*%v,N_Str)}))
    
    
  }
  
  return(list(ALL_COEFS=ALL_COEFS, BIC=BIC))	
}






mise_forme_coef_hiers_R = function(coefs, nstr)
{
  p        = length(coefs)/(nstr+1) 
  Coefbar  = coefs[1:p]
  Coeforme = Coefbar
  coefsred = coefs[-(1:p)]
  for (c in 1:nstr)
  {
    Coeforme = c(Coeforme, Coefbar + coefsred[((c-1)*p+1):(c*p)])
  }
  return(Coeforme)	
}




mise_forme_coef_hiers = function(coefs, nstr)
{
  p        = length(coefs)/(nstr+1) 
  Coefbar  = coefs[1:p]
  Coeforme = NULL
  coefsred = coefs[-(1:p)]
  for (c in 1:nstr)
  {
    Coeforme = c(Coeforme, Coefbar + coefsred[((c-1)*p+1):(c*p)])
  }
  return(Coeforme)	
}

mise_forme_finale<-function(vect,nstr)
{
  matrix(vect,nrow=nstr,byrow=T)
}



Calise_ITCPT = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK)
  nstr = length(ncs)
  n    = nrow(Feat)
  p    = ncol(Feat)
  Fin  = matrix(0, n, (p+1)*nstr)
  i0   = 1
  for (c in 1:nstr)
  {
    Fin[i0:(i0-1+ncs[c]), ((c-1)*(p+1)+1):(c*(p+1))] = cbind(rep(1, ncs[c]), Feat[i0:(i0-1+ncs[c]),])
    i0 										 = i0 + ncs[c]
  }
  return(Fin)
}


Hierise_ITCPT = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK) puis 1{C=1},...,1{C=K} puis X, puis 1
  Fin = cbind(Feat, Calise_ITCPT(Feat,ncs))
  return(Fin)
}




SepLogit_DataShared <-function(data,adap,standardize=FALSE){
  
  data     <-traitdata(data)
  p        <-ncol(data)
  ntot     <-nrow(data)
  N_Str    <-max(data[,p])
  X        <-data[,-p]
  X        <-as.matrix(X)
  Z        <-data[,p]
  p        <-ncol(X)
  ncs      <-table(Z)
  strata   <-1:N_Str

  
  Theta_BIC  = vector("list")
  
  for (kp in 1:p)
  { 

    Theta_BIC[[kp]]  = matrix(Inf, N_Str, p)
    
    Covars   = X[,-kp]
    p1       = ncol(Covars)
    respon   = X[,kp]
    
    if(standardize==TRUE){Covars=Covars%*%diag(1/((apply(Covars, 2, var))^0.5))}
    
    calX = Calise_ITCPT(Covars, ncs)
    
    
    if(adap==TRUE){calCoef  = glmnet(calX,respon,intercept=FALSE,family="binomial",alpha=1,lambda=0)$beta
    
    #if(adap==TRUE){calCoef  =glm(respon~calX-1, family="binomial")$coef
    
    tempo   = compute.alpha.beta.bar.DELTA(calCoef, p, N_Str, func=median_VV)
    
    pdsLam1 = 1/abs(tempo[[1]][-1])
    pdsLam2 = 1/abs(as.numeric(t(tempo[[2]])))
    pdsLam2[pdsLam2 == Inf] = 1e8}else{pdsLam1=pdsLam2=NULL}
    
    
    basic <- screen.DataSharedLasso.GLMs.ITCPT(Covars, respon, Z, p1, N_Str, family="binomial", Balanced= F, wObs= rep(1, length(respon)), ratio = 0.001, Nlam = 50, wLambda1=pdsLam1, wLambda2 = pdsLam2)

    bic = basic[[1]][,which.min(basic[[2]])]
    Theta_BIC[[kp]][,-kp] = mise_forme_finale(bic, N_Str)[,-1] 
    Theta_BIC[[kp]][,kp]  = mise_forme_finale(bic, N_Str)[,1] 


  }
  
  
  
  res <- list()
  res$MIN <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_BIC,p,N_Str)[[1]], N_Str, p)
  res$MAX <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_BIC,p,N_Str)[[2]], N_Str, p)

  return(res)
  }




SepLogit_Ref <-function(data,adap,standardize=FALSE)
  
{
  data     <-traitdata(data)
  p        <-ncol(data)
  ntot     <-nrow(data)
  N_Str    <-max(data[,p])
  X        <-data[,-p]
  X        <-as.matrix(X)
  Z        <-data[,p]
  p        <-ncol(X)
  ncs      <-table(Z)
  strata   <-1:N_Str

  
  Theta_BIC = vector("list")
  
  for (kp in 1:p)
  { 

    Theta_BIC[[kp]]  = matrix(Inf, N_Str, p)
    
    Covars   = X[,-kp]
    p1       = ncol(Covars)
    respon   = X[,kp]
    
    if(standardize==TRUE){Covars=Covars%*%diag(1/((apply(Covars, 2, var))^0.5))} 
    
    calX = Calise_ITCPT(Covars, ncs)
    calX =calX[,-(1:p)]
    
    #if(adap==TRUE){calCoef  = glmnet(calX,respon,intercept=FALSE,family="binomial",alpha=0,lambda=lam)$beta
    if(adap==TRUE){
      calCoef  =glm(respon~calX-1, family="binomial")$coef
      
      tempo   = compute.alpha.beta.bar.DELTA(c(rep(0,p),calCoef), p, N_Str, func=median_VV)
      
      pdsLam1 = 1/abs(tempo[[1]][-1])
      pdsLam2 = 1/abs(as.numeric(t(tempo[[2]])))
      pdsLam2[pdsLam2 == Inf] = 1e8}else{pdsLam1=pdsLam2=NULL}
    
    
    basic <- screen.RefLasso.GLMs.ITCPT(Covars, respon, Z, p1, N_Str, family="binomial", Balanced= F, wObs= rep(1, length(respon)), ratio = 0.001, Nlam = 50, wLambda1=pdsLam1, wLambda2 = pdsLam2)

    
    bic = basic[[1]][,which.min(basic[[2]])]
    Theta_BIC[[kp]][,-kp] = mise_forme_finale(bic, N_Str)[,-1] 
    Theta_BIC[[kp]][,kp]  = mise_forme_finale(bic, N_Str)[,1] 
    
  }
  
  
  
  res <- list()
  res$MIN <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_BIC,p,N_Str)[[1]], N_Str, p)
  res$MAX <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_BIC,p,N_Str)[[2]], N_Str, p)
  return(res)
}



median_VV = function(v)
{
  if (length(v)%%2 == 1){median(v)}
  else (median(c(-Inf, v)))
}




compute.alpha.beta.bar.DELTA = function(coefs.hier, p, N_Str, func=Mode)
{
  Mat.coef        = mise_forme_finale(coefs.hier, N_Str)
  alpha.beta.bar  = apply(Mat.coef[,1:p],2,func)
  DELTA           = Mat.coef[,1:p] - matrix(rep(alpha.beta.bar, N_Str),nrow=N_Str, byrow=T)
  return(list(alpha.beta.bar, DELTA))
}
