# Function qui construit les modèles fused pour des séquences de lam1 et lam2 "données"
# -> permet de screener parmi l'ensemble des modèles possibles et d'en renvoyer
# au maximum N1*N2 (si Nj représente la longueur de la séquence pour lamj)
transf_mat_sparse_format<-function(M,n)
{
	ind_nonzero<-which(M!=0)
	val_nonzero<-M[ind_nonzero]
	ii<-1+(ind_nonzero-1)%%nrow(M)
	jj<-1+(ind_nonzero-1)%/%nrow(M)
	Feat<-sparseMatrix(i=ii, j=jj, x=val_nonzero)
    val=Feat@x
    idx_i=Feat@i
    idx_p=Feat@p
    #val=ifelse(val<=1,val,1)
    Msparse=sparseMatrix(i=idx_i+1, p=idx_p, x=val, dims=c(n,nnn=ncol(Feat)))
	return(Msparse)
}

library(MatrixModels)
library(glmnet)


#### Avec le calcul des lambda2.max "lambda1-dépendants"

compute.lam2.max.noInt<-function(X,y,l1,pF,family = "gaussian" ,N_Str,adapt_weight_lam1=NULL,maxLam2,graph,wObs=NULL)
{
	
	pFused = pF 
	p      = pF/N_Str
	A = matrix(0,pFused,pFused)
	for (k in 1:pF)
	{
		Edges = graph$graphConn[[k]]
		A[k,Edges] = 1
	}
	diag(A) = rep(0,pF)	
	l2     = maxLam2

    b      = fusedlasso(X,y,lambda1=l1,lambda2=l2, family = family, wLambda1 =adapt_weight_lam1, graph=graph, verbose=FALSE,wObs=wObs,addIntercept=F)$beta ##ici
	nbmax  = sum(A)/2
    indeg = which(A==1, arr.ind=TRUE)
    forh   = sapply(1:nrow(indeg), FUN=function(iii){1*(b[indeg[iii,1]]==b[indeg[iii,2]])})
    h      = sum(unlist(forh))/2
    jjj    =0     
	while ( (h==nbmax) & (sum(b)!=0))
	{
		l2   = l2/2
		jjj  = jjj+1
		b    = fusedlasso(X,y,lambda1=l1,lambda2=l2, family = family, wLambda1 =adapt_weight_lam1, graph=graph, verbose=FALSE,wObs=wObs,addIntercept=F)$beta ##ici
        forh   = sapply(1:nrow(indeg), FUN=function(iii){1*(b[indeg[iii,1]]==b[indeg[iii,2]])})
    	h      = sum(unlist(forh))/2
	}
	if (jjj>0) {l2=l2*2}
	return(max(l2,1e-5))
}


constr_fused_Sparse_noInt<-function(X, y, pF, N_Str,family = "gaussian", adapt_weight_lam1=NULL, graph=NULL, ratio = 0.001, NN = 50, wObs = NULL, Threshold_LamMax1=0.1)
{
	# First, we compute lambda1 and lambda2 sequences
'	maxLams      = fusedlassoMaxLambdas(X, y, family = family, wLambda1 = adapt_weight_lam1, graph=graph, wObs = wObs,addIntercept=F)      ##ici
	maxLam1      = min(maxLams$maxLambda1,Threshold_LamMax1) 
	maxLam2      = maxLams$maxLambda2
	seq.lam1     = make_lambda_seq(maxLam1, ratio,NN)
	seq.lam2.max = sapply(seq.lam1,FUN = function(l1){compute.lam2.max.noInt(X,y,l1,pF,family = family,N_Str,adapt_weight_lam1,maxLam2,graph,wObs=NULL)})
	'
	
	S=adapt_weight_lam1
	if(is.null(S)){S=rep(1,nrow(X))}
	GL = glmnet(X,y, intercept=F, standardize = FALSE, family=family, penalty.factor = S) # rep(1,160))#
	
	if(is.null(S)){maxLam1 = max(GL$lambda)}else{
	  maxLam1 = max(GL$lambda)*(nrow(X)/sum(S))}
	maxLam2 =1000
	seq.lam1    = make_lambda_seq(maxLam1, ratio,NN)
	seq.lam2.max = sapply(seq.lam1,FUN = function(l1){compute.lam2.max.noInt(X,y,l1,pF,family = family,N_Str,adapt_weight_lam1,maxLam2,graph,wObs=NULL)})
	
	
	
	Betas<-vector("list", NN)
	Inters<-vector("list", NN)
	
	betaStart1<-NULL; interceptStart1<-NULL
	betaStart2<-NULL; interceptStart2<-NULL
	klam2<-0; klam1<-0
	
	#######################################
	#   on applique le code d'Holger, pour
	# chaque valeur du couple lambda_1,lambda_2
	# (pour gagner du temps, on commence par les 
	# valeurs élevées de lambda_2 et/ou lambda_1
	# et on utilise les solutions obtenues comme valeurs
	# initiales lorsqu'on diminue les
	# valeurs de lambda_1 et/ou lambda_2) 	     
	#######################################

	
	for (lam1 in seq.lam1)
	{
			klam2   = 0
			klam1   = klam1+1
			seqlam2 = make_lambda_seq(seq.lam2.max[klam1], ratio,NN)
			RES<-fusedlasso(X, y, lambda2=seqlam2[1],lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, 
					graph = graph, betaStart = betaStart2, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
			betaStart2<-as.numeric(RES$beta); 
			betaStart1<-betaStart2;
			TEMP<-as.numeric(RES$beta)

			for (lam2 in seqlam2[-1])
			{
				klam2<-klam2+1
				RES<-fusedlasso(X, y, lambda2=lam2,lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, 
				graph = graph, betaStart = betaStart1, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
				betaStart1<-as.numeric(RES$beta); 
				TEMP<-cbind(TEMP,as.numeric(RES$beta))
			}
			Betas[[klam1]]<-TEMP	
			rm(TEMP)
	}	
	return(list(Betas))	
}



compute.beta.bar.DELTA.Fused = function(coefs.fused, p, N_Str)
{
	Mat.coef        = mise_forme_finale(coefs.fused, N_Str)
	beta.bar        = apply(Mat.coef[,1:p],2,Mode)
	DELTA           = Mat.coef[,1:p] - matrix(rep(beta.bar, N_Str),nrow=N_Str, byrow=T)
	return(list(beta.bar, DELTA))
}
		
sel.BIC.relaxed.Fused = function(ALL_COEFS, X, y, Z, p, N_Str, family = "gaussian", pen.ridge = 1e-7, pen=log(length(y)))
{		

	comp.BIC  = function(coefs.hier) 
	{
		DECOMP = compute.beta.bar.DELTA.Fused(coefs.hier,p,N_Str)
		comp.BIC.relaxed.logisticlinear(DECOMP[[1]],DECOMP[[2]], X, y, Z, p, N_Str, family = family, pen.ridge = pen.ridge, pen=log(length(y)))
	} 
	
	TEST_ALL  = apply(ALL_COEFS, 2,comp.BIC)
	BICs      = unlist(lapply(TEST_ALL, FUN=function(x){x$BICc}))
	
	return(mise_forme_coef_hiers(TEST_ALL[[which.min(BICs)]]$coef.relax,N_Str))
	
}


K_Fold_constr_fused_Sparse_noInt<-function(X, y, pF, N_Str, family = "gaussian", Kfold=5, adapt_weight_lam1=NULL, graph=NULL, ratio = 0.001, NN = 50, wObs = NULL, Threshold_LamMax1=0.1)
{
	# First, we compute lambda1 and lambda2 sequences
	maxLams      = fusedlassoMaxLambdas(X, y, family = family, wLambda1 = adapt_weight_lam1, graph=graph, wObs = wObs,addIntercept=F)      ##ici
	maxLam1      = min(maxLams$maxLambda1,Threshold_LamMax1) 
	maxLam2      = maxLams$maxLambda2
	seq.lam1     = make_lambda_seq(maxLam1, ratio,NN)
	seq.lam2.max = sapply(seq.lam1,FUN = function(l1){compute.lam2.max.noInt(X,y,l1,pF,N_Str,adapt_weight_lam1,maxLam2,graph,wObs=NULL)})
	
	
	Betas<-vector("list", NN)
	Inters<-vector("list", NN)
	
	betaStart1<-NULL; interceptStart1<-NULL
	betaStart2<-NULL; interceptStart2<-NULL
	klam2<-0; klam1<-0
	
	#######################################
	#   on applique le code d'Holger, pour
	# chaque valeur du couple lambda_1,lambda_2
	# (pour gagner du temps, on commence par les 
	# valeurs élevées de lambda_2 et/ou lambda_1
	# et on utilise les solutions obtenues comme valeurs
	# initiales lorsqu'on diminue les
	# valeurs de lambda_1 et/ou lambda_2) 	     
	#######################################
	Betas = NULL
	for (lam1 in seq.lam1)
	{
			klam2   = 0
			klam1   = klam1+1
			seqlam2 = make_lambda_seq(seq.lam2.max[klam1], ratio,NN)
			RES<-fusedlasso(X, y, lambda2=seqlam2[1],lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, 
					graph = graph, betaStart = betaStart2, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
			betaStart2<-as.numeric(RES$beta); 
			betaStart1<-betaStart2;
			TEMP<-as.numeric(RES$beta)

			for (lam2 in seqlam2[-1])
			{
				klam2<-klam2+1
				RES<-fusedlasso(X, y, lambda2=lam2,lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, 
				graph = graph, betaStart = betaStart1, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
				betaStart1<-as.numeric(RES$beta); 
				TEMP<-cbind(TEMP,as.numeric(RES$beta))
			}
			Betas = cbind(Betas,TEMP)	
			rm(TEMP)
	}	
	
	ncv         = length(y)/N_Str
	RndOBS0     = sample(1:ncv,ncv)
	err_pred_cv = 0
	for (kk in 1:Kfold)
	{		
		All_Coef_tr = NULL
		
		RndOBS      = outer(RndOBS0[((kk-1)*floor(ncv/Kfold)+1):(kk*floor(ncv/Kfold))], (0:(N_Str-1))*ncv, "+")	
		
		Part_test   = sort(as.numeric(RndOBS))
		
		Betas_tr = NULL
		
		klam2<-0; klam1<-0
		for (lam1 in seq.lam1)
		{
			klam2   = 0
			klam1   = klam1+1
			seqlam2 = make_lambda_seq(seq.lam2.max[klam1], ratio,NN)
			RES<-fusedlasso(X[-Part_test, ], y[-Part_test], lambda2=seqlam2[1],lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, graph = graph, betaStart = betaStart2, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
			betaStart2<-as.numeric(RES$beta); 
			betaStart1<-betaStart2;
			TEMP<-as.numeric(RES$bketa)

			for (lam2 in seqlam2[-1])
			{
				klam2<-klam2+1
				RES<-fusedlasso(X[-Part_test, ], y[-Part_test], lambda2=lam2,lambda1=lam1, family = family, wLambda1 = adapt_weight_lam1, graph = graph, betaStart = betaStart1, addIntercept=F, verbose=FALSE,wObs=wObs) ##ici
				betaStart1<-as.numeric(RES$beta); 
				TEMP<-cbind(TEMP,as.numeric(RES$beta))
			}
			Betas_tr = cbind(Betas_tr,TEMP)	
			rm(TEMP)
		}
		for_errtest = (y[Part_test]- X[Part_test,]%*%Betas_tr)^2		
		err_pred_cv = err_pred_cv + colSums(for_errtest)		
	}
	
	return(Betas[, which.min(err_pred_cv)])	
}


cree_graph_compute_weights_forVOs = function(calCoef, p, N_Str)
{
	which_median_VV = function(v)
	{
		which(v==median_VV(v))
	}
	pFused 			   = p*N_Str
	graphConn.star     = vector("list", pFused)
	graphWeight.star   = vector("list", pFused)
	Inds_stratcomp_ref = NULL
	for (k in 1:p)
	{
		# On cherche la strate qui correspond à la strate de référence, pour chaque feature
		# sert pour le star-graph (VO1 & VO2)
		# mais aussi pour la pénalité Lasso (VO1 : seule celle-ci est pénalisée)	
		ind_coef_comp_k    = (0:(N_Str-1))*p  + k	
		ind_stratcomp_ref  = ind_coef_comp_k[which_median_VV(calCoef[ind_coef_comp_k])]
		Inds_stratcomp_ref = c(Inds_stratcomp_ref,ind_stratcomp_ref)
	
		for (jk in ind_coef_comp_k)
		{
			if (jk == ind_stratcomp_ref)
			{
				Edges = ind_coef_comp_k[-which(ind_coef_comp_k==ind_stratcomp_ref)] 

			}else
			{
				Edges = ind_stratcomp_ref
			}
			graphConn.star[[jk]]   = as.integer(Edges)
			graphWeight.star[[jk]] = 1/abs(calCoef[jk]-calCoef[Edges])
		}
	}				

	graph.star       = list(graphConn=graphConn.star, graphWeight=graphWeight.star)

# 2 vecteurs de poids, pour VO1 et VO2
	wLambda1_VO1 = rep(0,pFused)
	wLambda1_VO1[Inds_stratcomp_ref] = 1/abs(calCoef[Inds_stratcomp_ref])
	wLambda1_VO2 = 1/abs(calCoef)
	
	return(list(graph.star=graph.star, wLambda1_VO1=wLambda1_VO1, wLambda1_VO2=wLambda1_VO2))
}


cree_graph_compute_weights_cliques = function(calCoef=NULL, p, N_Str)
{
	pFused = p*N_Str
	graphConn <- vector("list", pFused )
	if (!is.null(calCoef)) {graphWeight.adapt <- vector("list", pFused)}
	graphWeight <- vector("list", pFused )
	for (k in 1:pFused)
	{
		kmodp            = ifelse(k%%p==0,p,k%%p)
		Edges            = kmodp+(0:(N_Str-1))*p
		graphConn[[k]]   = as.integer(Edges[-(which(Edges==k))])
		graphWeight[[k]] = rep(1,length(Edges)-1)
		if (!is.null(calCoef)) {graphWeight.adapt[[k]]<-1/abs(calCoef[k]-calCoef[Edges[-(which(Edges==k))]])}
	}				

	## 2 graphs : with & without adaptive weights
	graph       = list(graphConn=graphConn, graphWeight=graphWeight)
	if (!is.null(calCoef)) 
	{
		graph.adapt = list(graphConn=graphConn, graphWeight.adapt=graphWeight.adapt)
		wLambda1 	= 1/abs(calCoef)
	}else
	{
		graph.adapt = NULL
		wLambda1 	= NULL
	}
	
	return(list(graph=graph, graph.adapt=graph.adapt,wLambda1=wLambda1))
}

