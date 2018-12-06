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



CPX <- function(x)
{
  length(unique(x[x!=0]))
}


summarize_forMatrices = function(forMatrices, p,N_Str)
{
  
  Matrix_f_min = Matrix_f_max = vector("list")
  for (j1 in 1:p)
  {
    Matrix_f_min[[j1]] = Matrix_f_max[[j1]] = matrix(NA, N_Str, p)
  }
  
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {  
      o=0
      for (k in 1:N_Str) {
        if(is.na(forMatrices[[j1]][k,j2])== T || is.na(forMatrices[[j2]][k,j1])==T){o=1}}
      
      if (j1 != j2)
      {
        if(o==0){
          Vect1 = forMatrices[[j1]][,j2]
          Vect2 = forMatrices[[j2]][,j1]
          
          Complex1 = CPX(Vect1)
          Complex2 = CPX(Vect2)
          
          if (Complex1 < Complex2)
          {
            Matrix_f_min[[j1]][,j2] = Matrix_f_min[[j2]][,j1] = Vect1
            Matrix_f_max[[j1]][,j2] = Matrix_f_max[[j2]][,j1] = Vect2
          }else
          {
            Matrix_f_min[[j1]][,j2] = Matrix_f_min[[j2]][,j1] = Vect2
            Matrix_f_max[[j1]][,j2] = Matrix_f_max[[j2]][,j1] = Vect1
          }
        }else{
          Matrix_f_min[[j1]][,j2] = Matrix_f_max[[j1]][,j2] = forMatrices[[j1]][,j2]
          Matrix_f_min[[j2]][,j1] = Matrix_f_max[[j2]][,j1] = forMatrices[[j2]][,j1]}
      }else
      {
        Matrix_f_max[[j1]][,j1]=Matrix_f_min[[j1]][,j1]=forMatrices[[j1]][,j1]
      }
    }
  }
  return(list(Matrix_f_min, Matrix_f_max))
}


MiseForme_Matrices_f_maxmin = function(Matrix_f, N_Str, p)
{
  Matrices_FIN = vector("list")
  for (knstr in 1:N_Str)
  {
    Matrices_FIN[[knstr]] = matrix(NA, p, p)
    for (kp in 1:p)
    {
      Matrices_FIN[[knstr]][,kp] = Matrix_f[[kp]][knstr,]
    }
  }
  Matrices_FIN<-array(unlist(Matrices_FIN),c(p,p,N_Str))
  return(Matrices_FIN)	
}
