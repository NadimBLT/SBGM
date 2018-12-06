make_lambda_seq<-function(lambda_max, ratio,N)
{
  lambda_seq<-exp(seq(log(lambda_max),log(lambda_max*ratio), length.out=N))		
}




Calise = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK)
  nstr = length(ncs)
  n    = nrow(Feat)
  p    = ncol(Feat)
  Fin  = matrix(0, n, p*nstr)
  i0   = 1
  for (c in 1:nstr)
  {
    Fin[i0:(i0-1+ncs[c]), ((c-1)*p+1):(c*p)] = Feat[i0:(i0-1+ncs[c]),]
    i0 										 = i0 + ncs[c]
  }
  return(Fin)
}

SepLogit_Fused <- function(data,adap,standardize=FALSE){
  
  data     <-traitdata(data)
  p        <-ncol(data)
  ntot     <-nrow(data)
  Nlam     <-50
  N_Str    <-max(data[,p])
  X        <-data[,-p]
  X        <-as.matrix(X)
  Z        <-data[,p]
  p        <-ncol(X)
  ncs      <-table(Z)
  strata   <-1:N_Str
  
  Theta_Fused =  vector("list")
  
  for (kp in 1:p)
  { 
    
    Theta_Fused[[kp]]= matrix(Inf, N_Str, p)
    
    Covars   = X[,-kp]
    respon   = X[,kp]
    
    if(standardize==TRUE){Covars=Covars%*%diag(1/((apply(Covars, 2, var))^0.5))}
    
    Covars   = cbind(rep(1,ntot),Covars)
    calX = Calise(Covars, ncs)
    
    if(adap==TRUE){calCoef  = glm(respon~calX-1, family="binomial")$coef}else{calCoef=NULL}
    
    CGCWC = cree_graph_compute_weights_cliques(calCoef=calCoef, p, N_Str)
    graph = CGCWC$graph
    pFused = p*N_Str
    if(adap==TRUE){
      SCREEN  = constr_fused_Sparse_noInt(calX,respon,pFused, N_Str,family="binomial", adapt_weight_lam1=1/abs(calCoef), graph= graph)}else{
        SCREEN  = constr_fused_Sparse_noInt(calX,respon,pFused, N_Str,family="binomial", adapt_weight_lam1=NULL, graph= graph)}
    
    
    kkk         = 0
    Coef.Fused  = matrix(NA, pFused, Nlam^2)
    
    for (kk1 in 1:Nlam){
      for (kk2 in 1:Nlam){
        kkk                    = kkk+1
        Coef.Fused[,kkk]       = SCREEN[[1]][[kk1]][,kk2]
      } 
    }
    
    BIC=rep(NA,Nlam^2)
    
    
    for(j in 1:Nlam^2){
      LL=0
      for (i in 1:length(respon)) {
        
        LL = LL + (respon[i]*(t(Coef.Fused[,j])%*%calX[i,]) - log(1+exp(t(Coef.Fused[,j])%*%calX[i,])))
      }
      BIC[j]= -2*LL + length(unique(as.vector(Coef.Fused[,j]))!=0)*log(length(respon))
    }
    
    
    Fused.bic <- Coef.Fused[,which.min(BIC)]
    
    Theta_Fused[[kp]][,-kp] = mise_forme_finale(Fused.bic, N_Str)[,-1]
    Theta_Fused[[kp]][,kp] = mise_forme_finale(Fused.bic, N_Str)[,1]
    
    
    
    
  }
  
  res <- list()
  res$MIN <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_Fused,p,N_Str)[[1]], N_Str, p)
  res$MAX <-MiseForme_Matrices_f_maxmin(summarize_forMatrices(Theta_Fused,p,N_Str)[[2]], N_Str, p)
  
  
  
  return(res)
}




