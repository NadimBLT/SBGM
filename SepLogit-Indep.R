

IndepLasso_ITCPT_BIC <-function(X, Y, Z, p=ncol(X), N_str=max(Z), family = "binomial", Nlam=50, ratio=0.001,W=FALSE){
  
  RES <- lapply(1:N_str, function(k){
    BIC=NULL
    ALL_Coefs=NULL
    
    Y_k <- Y[which(Z==k)]
    X_k <- X[which(Z==k),]

    wObs= rep(1, length(Y_k))
      
    if (W==FALSE) {
      W=rep(1,p)
    }else{
      W =glm(Y_k~X_k-1, family="binomial")$coef
    }
    
    
    MatDiagWeightsInv = sparseMatrix(i = 1:(p), j = 1:(p) ,x= 1/W)
    X_k=X_k%*%MatDiagWeightsInv 
    
    
    
    mod0  <- glmnet(X_k, Y_k, weights= wObs, standardize=F, family=family, nlambda=Nlam, lambda.min.ratio=ratio)
      
      dev0         = mod0$nulldev*(1-mod0$dev.ratio)
      Devf         = dev0 - mod0$nulldev
      BIC          = Devf + mod0$df*log(length(Y_k))

      #ALL_Coefs   = mod0$beta
      bbeta	  = rbind(mod0$a0, as.matrix(mod0$beta))
      bbeta		= as.matrix(bbeta)
      
      MatDiagWeightsInv_ITCPT = Matrix(cbind(c(1, rep(0, nrow(MatDiagWeightsInv))), rbind(rep(0, ncol(MatDiagWeightsInv)), as.matrix(MatDiagWeightsInv))), sparse=T)
      ALL_Coefs = cbind(ALL_Coefs, apply(bbeta,2,function(v){as.vector(MatDiagWeightsInv_ITCPT%*%v)}))
      
      
      return(list(ALL_Coefs,
                  BIC))
    }
  )
  return(RES)
}








SepLogit_Indep <-function(data,standardize=FALSE,adap=FALSE)
  
{
  data     <-traitdata(data)
  p        <-ncol(data)
  ntot     <-nrow(data)
  N_str    <-max(data[,p])
  X        <-data[,-p]
  X        <-as.matrix(X)
  Z        <-data[,p]
  p        <-ncol(X)
  ncs      <-table(Z)
  strata   <-1:N_str
  
  Theta_IndepLasso_BIC  = vector("list")
  
  for (j in 1:p)
  { 
    
    Theta_IndepLasso_BIC[[j]] =matrix(NA, N_str, p)
    
    Covars   = X[,-j]
    respon   = X[,j]
    
    if(standardize==TRUE){Covars=Covars%*%diag(1/((apply(Covars, 2, var))^0.5))}
    
    
    RES=IndepLasso_ITCPT_BIC(X=Covars, Y=respon, Z, p=ncol(Covars), N_str=N_str, family = "binomial", Nlam=50, ratio=0.001,W=adap)
    
    bic2=sapply(1:N_str, function(k){return(RES[[k]][[1]][,which.min(RES[[k]][[2]])])})
    Theta_IndepLasso_BIC[[j]][,-j] = t(bic2)[,-1] 
    Theta_IndepLasso_BIC[[j]][,j]  = t(bic2)[,1] 
    
    
  }
  res <- list()
  res$MIN <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_IndepLasso_BIC,p,N_str)[[1]], N_str, p)
  res$MAX <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_IndepLasso_BIC,p,N_str)[[2]], N_str, p) 
  
  return(res) 
}











