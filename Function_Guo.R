


Function_Guo    = function(simData_K)
{
  
	non_zero     <-recup_nonzero(simData_K)
	lambda.max   <-choix_lambda_max(simData_K)
	lambda.min   <-0.001
	rhoVec       <-exp(seq(log(lambda.max), log((lambda.min*lambda.max)), length.out = 50))
	Theta_init   <-calc_Theta(simData_K,rhoVec)        
	Theta_in     <-Theta_init
	p            <-dim(Theta_init)[1]
	k            <-dim(Theta_init)[3]
	nbe.lam      <-length(rhoVec)
	Theta_update <-array(NA,c(p,p,k,nbe.lam))
	rhoVec_update<-array(NA,c(p,p,nbe.lam))

	for (j in 1:nbe.lam)                   
	{
	    pseudo1_n<-c()
  		pseudo2_n<-c()
  		max_vec  <-c()
  		conv<-FALSE
  		maxItera<-0
  		while (conv==FALSE)
  		{
    		maxItera          <-maxItera+1
    		afac              <-rhoVec[j]/threshold(sqrt(apply(abs(Theta_init[,,,j]),c(1,2),sum)))   # function threshold 
    		diag(afac)        <-0
    		rhoVec_update[,,j]<-afac
    		for (i in 1:k)  
    		{
      			Theta_update[,,i,j]<-BMNPseudo.single(choix_groupe(simData_K,i),rhoVec_update[,,j],ThetaStart=Theta_init[,,i,j])$Theta 
    		}    
    		if ((max(abs(Theta_init[,,,j]-Theta_update[,,,j]))<1e-5)==TRUE | maxItera==100)
    		{
      			conv<-TRUE
      			print(maxItera)
    		}
    		Theta_init[,,,j]<-Theta_update[,,,j]      
  		}
	} 

	n             <-nrow(simData_K)
	Theta_update  <-traitTheta(Theta_update,non_zero,simData_K)
	p             <-dim(Theta_update)[1]
	Theta_nonbiais<-array(NA,c(p,p,k,nbe.lam))
	err_estim_upd <-c(rep(0,nbe.lam))
	acc_sup       <-c(rep(0,nbe.lam))
	err_estim_nb  <-c(rep(0,nbe.lam))
	bic           <-c(rep(0,nbe.lam))
	bicVV         <-c(rep(0,nbe.lam))
	aic           <-c(rep(0,nbe.lam))
	aicVV         <-c(rep(0,nbe.lam))
	a1            <-c()
	a2            <-c()
	a3            <-c()

	for (j in 1:nbe.lam)
	{
  		sum_pL    <-0
  		coef_nonul<-0
  		for (i in 1:k)
  		{
		    conc.constr          <-ifelse(Theta_update[,,i,j]==0,Inf,0)
    		Theta_nonbiais[,,i,j]<-BMNPseudo.single(choix_groupe(simData_K,i),conc.constr,ThetaStart=Theta_update[,,i,j])$Theta
		    coef_nonul           <-coef_nonul+sum(Theta_nonbiais[,,i,j][upper.tri(Theta_nonbiais[,,i,j])]!=0)
		    
		    dota                 <-choix_groupe(simData_K,i)   
    		for (l in 1:nrow(dota))
    		{
      			for (m in 1:ncol(dota))
      			{
        			z        <-2*dota[l,m]-1
			        dota1    <-dota
			        dota1[,m]<-1
        			sum_pL   <-sum_pL-log(1+exp(-z*dota1[l,]%*%Theta_nonbiais[m,,i,j]))
      			}
    		}    		
    	}    	
  	
  	
  	bic[j]  <-(-2)*sum_pL+log(n)*coef_nonul
  	bicVV[j]<-(-2)*sum_pL/2+log(n)*coef_nonul
  	aic[j]  <-(-2)*sum_pL+2*coef_nonul
  	aicVV[j]<-(-2)*sum_pL/2+2*coef_nonul
}
	#return(list(BIC=Theta_nonbiais[,,,which.min(bic)], Theta_nonbiais[,,,which.min(bicVV)], Theta_nonbiais, bic, bicVV))  
	return(list(BIC=Theta_nonbiais[,,,which.min(bic)]))
}



