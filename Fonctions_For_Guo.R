traitdata<-function(data) {       # ordonne en fonction du groupe et enleve les colonnes "nulles"
  p    <-ncol(data)
  data <-data[order(data[,p]),]
  duta <-data[,-p]
  categ<-data[,p]
  if (min(duta)==0) {duta<-2*as.matrix(duta)-1}
  
  nonzero=unique(which(colSums(duta==-1)>0 & colSums(duta==1)>0))
  dota=duta[,nonzero]
  
  if (min(dota)==-1) {dota<-(dota+1)/2}
  
  data<-cbind(dota,k=categ)
  return(data)
}

traitTheta<-function(Theta,nonzero,data)
{
  p<-ncol(data)-1
  k<-dim(Theta)[3]
  nbe_lam<-dim(Theta)[4]
  Theta_f<-array(0,c(p,p,k,nbe_lam))
  Theta_f[nonzero,nonzero,,]<-Theta
  return(Theta_f)
}

recup_nonzero<-function(data) {       
  p    <-ncol(data)
  data <-data[order(data[,p]),]
  duta <-data[,-p]
  duta<-2*as.matrix(duta)-1
  nonzero=unique(which(colSums(duta==-1)>0 & colSums(duta==1)>0))
  return(nonzero)
}

choix_groupe<-function(data,i){  
  if (max(data[,ncol(data)])<i)
  {
    stop("wrong category argument")
  }
  data<-traitdata(data)
  p   <-ncol(data)
  dota<-data[,-p]
  dota<-as.matrix(dota)
  
  dota<-dota[which(data[,p]==i),]
  return(dota)
}

choix_lambda_max<-function(data)
{     
  k         <-max(data[,-ncol(data)])
  lambda.max<-0
  for (i in 1:k)
  {
    dota<-choix_groupe(data,i)
    n   <-nrow(dota)
    dota<-2*as.matrix(dota)-1
    p   <-ncol(dota)
    
    for (pvar1 in 1:p)
    {
      B         <-matrix(dota[,pvar1],n,1)
      Bp1       <-matrix(rep(B,p-1),n,p-1)
      A         <-Bp1*dota[,-pvar1]
      mp        <-sum(B==1)
      mm        <-sum(B==-1)
      lambdamax1<-max(t(A)%*%(1/(1+exp(log(mp/mm)*B))))/n
      lambda.max<-max(lambdamax1,lambda.max)
    }
  }  
  return(lambda.max) 
}

threshold<-function(theta_k)
{
  theta_k[theta_k<1e-10] <-0
  return(theta_k)
}  

calc_Theta<-function(data,rhoVec)                                                  #verbose=T pour afficher lambda
{
  data      <-traitdata(data)
  p         <-ncol(data)
  dota      <-data[,-p]
  dota      <-as.matrix(dota)
  nbe.lam   <-length(rhoVec)
  k         <-max(data[,p])
  pseudoPath<-array(NA,c(p-1,p-1,k,nbe.lam))
  
  for (i in 1:k)
  {
    pseudoPath[,,i,]<-array(unlist(BMNPseudo(choix_groupe(data,i),rhoVec,ThetaStart=NULL)$ThetaList),c(p-1,p-1,1,nbe.lam))
  }
  
  return(pseudoPath)
}
