library(MASS)
library(Matrix)
#library(pheatmap)
library(regress)

MCMC_heritability<-function(phenotype,G,num_samples=50000){

    # 10-16-2019: To make prior on heritability to be uniform distribution on
    #             [0,1], I change a from -1 to 1 and b from 0 to 0.001

    # 10-17-2019: Changing a from -1 to 1 on 10-16-2019 gave very poor
    #             estimates of heritability and highly autocorrelated
    #             samples in MCMC chain, so I reverted to a = -1.
    #             The estimates seemed to improve.

    a<- -1
    b<-0.001

    y<-phenotype
    n<-length(y)
 
    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    X<-t(chol(G))

    XtX<-t(X)%*%X
    In<-diag(1,n)

    sigma_sq_gamma<-rep(NA,num_samples)
    sigma_sq<-rep(NA,num_samples)
    h_sq<-rep(NA,num_samples)
    mu<-rep(NA,num_samples)
    beta_mat<-matrix(NA,n,num_samples)

    sigma_sq_gamma[1]<-1/2*var(y)
    sigma_sq[1]<-1/2*var(y)
    numer<-sum(diag(sigma_sq_gamma[1]*G))
    denom<-numer+n*sigma_sq[1]
    h_sq[1]<-numer/denom
    mu[1]<-mean(y)
    beta_mat[,1]<-colSums(t(solve(XtX+In)%*%t(X))*(y-mu[1]))

    for(i in 2:num_samples){
        tau_gamma<-rgamma(1,shape=n/2+a,rate=1/2*sum(beta_mat[,i-1]^2)+b)
        sigma_sq_gamma[i]<-1/tau_gamma

        tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y-mu[i-1]-colSums(t(X)*beta_mat[,i-1]))^2)+b)
        sigma_sq[i]<-1/tau

        numer<-sum(diag(sigma_sq_gamma[i]*G))
        denom<-numer+n*sigma_sq[i]
        h_sq[i]<-numer/denom

        mu[i]<-rnorm(1,mean(y-colSums(t(X)*beta_mat[,i-1])),sqrt(sigma_sq[i]/n))

        beta_cov<-sigma_sq[i]*solve(XtX+sigma_sq[i]/sigma_sq_gamma[i]*In)
        beta_mean<-colSums(t((1/sigma_sq[i])*beta_cov%*%t(X))*(y-mu[i]))

        beta_mat[,i]<-mvrnorm(1,beta_mean,beta_cov)
        if((num_samples-i)%%1000==0){
          cat(i,"/", num_samples, "samples.\n")
        }
    }
    return(list(beta=beta_mat,mu=mu,sigma_sq=sigma_sq,sigma_sq_gamma=sigma_sq_gamma,h_sq=h_sq))
}


MCMC_heritability_regression<-function(form,dat,G,num_samples=50000){

    # 10-16-2019: To make prior on heritability to be uniform distribution on
    #             [0,1], I change a from -1 to 1 and b from 0 to 0.001

    # 10-17-2019: Changing a from -1 to 1 on 10-16-2019 gave very poor
    #             estimates of heritability and highly autocorrelated
    #             samples in MCMC chain, so I reverted to a = -1.
    #             The estimates seemed to improve.

    a<- -1
    b<-0.001

    fit_lm<-lm(form,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm)))
    Z<-model.matrix(fit_lm)
    cov_names<-dimnames(model.matrix(fit_lm))[[2]]

    ZtZinv<-solve(t(Z)%*%Z)
    ZtZinvZt<-ZtZinv%*%t(Z)

    p<-dim(Z)[2]
    n<-length(y)

    if(!is.null(fit_lm$na.action)){
      G<-G[-fit_lm$na.action,-fit_lm$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    X<-Matrix(t(chol(G)),sparse=TRUE)

    XtX<-Matrix(t(X)%*%X,sparse=TRUE)
    In<-Matrix(diag(1,n),sparse=TRUE)

    sigma_sq_gamma<-rep(NA,num_samples)
    sigma_sq<-rep(NA,num_samples)
    h_sq<-rep(NA,num_samples)
    alpha_mat<-matrix(NA,p,num_samples)
    beta_mat<-matrix(NA,n,num_samples)

    sigma_sq_gamma[1]<-1/2*var(fit_lm$resid)
    sigma_sq[1]<-1/2*var(fit_lm$resid)
    numer<-sum(diag(sigma_sq_gamma[1]*G))
    denom<-numer+n*sigma_sq[1]
    h_sq[1]<-numer/denom
    alpha_mat[,1]<-colSums(t(ZtZinvZt)*y)
    beta_mat[,1]<-colSums(t(chol2inv(chol(XtX+In))%*%t(X))*(y-colSums(t(Z)*alpha_mat[,1])))


    #sigma_sq_gamma<-NULL
    #sigma_sq<-NULL
    #h_sq<-NULL
    #alpha_mat<-NULL
    #beta_mat<-NULL

    #sigma_sq_gamma<-c(sigma_sq_gamma,1/2*var(fit_lm$resid))
    #sigma_sq<-c(sigma_sq,1/2*var(fit_lm$resid))
    #numer<-sum(diag(sigma_sq_gamma[1]*G))
    #denom<-numer+n*sigma_sq[1]
    #h_sq<-c(h_sq,numer/denom)
    #alpha_mat<-cbind(alpha_mat,colSums(t(ZtZinvZt)*y))
    #beta_mat<-cbind(beta_mat,colSums(t(chol2inv(chol(XtX+In))%*%t(X))*(y-colSums(t(Z)*alpha_mat[,1]))))


    for(i in 2:num_samples){
        tau_gamma<-rgamma(1,shape=n/2+a,rate=1/2*sum(beta_mat[,i-1]^2)+b)
        sigma_sq_gamma[i]<-1/tau_gamma
        #sigma_sq_gamma<-c(sigma_sq_gamma,1/tau_gamma)


        tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y-colSums(t(Z)*alpha_mat[,i-1])-colSums(t(X)*beta_mat[,i-1]))^2)+b)
        sigma_sq[i]<-1/tau
        #sigma_sq<-c(sigma_sq,1/tau)

        numer<-sum(diag(sigma_sq_gamma[i]*G))
        denom<-numer+n*sigma_sq[i]
        h_sq[i]<-numer/denom
        #h_sq<-c(h_sq,numer/denom)

        alpha_cov<-sigma_sq[i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y-colSums(t(X)*beta_mat[,i-1])))

        alpha_mat[,i]<-mvrnorm(1,alpha_mean,alpha_cov)
        #alpha_mat<-cbind(alpha_mat,mvrnorm(1,alpha_mean,alpha_cov))

        beta_cov<-sigma_sq[i]*chol2inv(chol(XtX+sigma_sq[i]/sigma_sq_gamma[i]*In))
        beta_mean<-colSums(t((1/sigma_sq[i])*beta_cov%*%t(X))*(y-colSums(t(Z)*alpha_mat[,i])))

        #beta_mat[,i]<-mvrnorm(1,beta_mean,beta_cov)
        beta_mat[,i]<-beta_mean+as.vector(t(chol(beta_cov))%*%rnorm(n))
        #beta_mat<-cbind(beta_mat,beta_mean+as.vector(t(chol(beta_cov))%*%rnorm(n)))

        if((num_samples-i)%%1000==0){cat(num_samples-i,' samples remaining.\n')}
        #print(i)
    }
    dimnames(alpha_mat)[[1]]<-cov_names
    return(list(beta=beta_mat,alpha=alpha_mat,sigma_sq=sigma_sq,sigma_sq_gamma=sigma_sq_gamma,h_sq=h_sq))
}



MCMC_heritability_mixed_genetic_and_non_genetic<-function(fixed_form,random_form,genetic_form,dat,G,num_chains=3,num_samples=10000,hyperpars_genetic=c(1,1),hyperpars_random=c(1,1),hyperpars_error=c(1,1)){

    # 12-18-2019: I set default prior on precision parameter to be
    #             Gamma(1,1). The implied prior on heritability is
    #             Uniform([0,1]).

    # 07-27-2021: For now, only one non-genetic random effect is allowed in random form.

    a_g<-hyperpars_genetic[1]
    b_g<-hyperpars_genetic[2]

    a_r<-hyperpars_random[1]
    b_r<-hyperpars_random[2]

    a_e<-hyperpars_error[1]
    b_e<-hyperpars_error[2]

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))
    form3<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(genetic_form)[2]))

    full_form<-as.formula(paste(paste(c(as.character(fixed_form)[2],as.character(fixed_form[1]),as.character(fixed_form)[3]),collapse=""),"+",as.character(random_form)[2],"+",as.character(genetic_form)[2],collapse=""))
    missing_vals<-lm(full_form,data=dat)$na.action

    if(!is.null(missing_vals)){
      dat <- dat[-missing_vals,]
    }
    
    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)
    fit_lm3<-lm(form3,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    scl_fctr<-1/sd(y)
    y_scaled<-y*scl_fctr
    Z<-model.matrix(fit_lm1)
    R<-model.matrix(fit_lm2)
    W_temp<-model.matrix(fit_lm3)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]
    cov_names3<-dimnames(model.matrix(fit_lm3))[[2]]
    

    ZtZinv<-chol2inv(chol(t(Z)%*%Z))
    ZtZinvZt<-ZtZinv%*%t(Z)
    RtR<-t(R)%*%R

    q<-length(W)
    p<-ncol(Z)
    m<-ncol(R)
    n<-length(y_scaled)

    if(!is.null(missing_vals)){
      G<-G[-missing_vals,-missing_vals]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    svdG<-svd(G)

    X<-Matrix(svdG$u%*%diag(sqrt(svdG$d)),sparse=TRUE)

    XtX<-Matrix(t(X)%*%X,sparse=TRUE)
    In<-Matrix(diag(1,n),sparse=TRUE)

    WX<-list()
    for(j in 1:q){
      WX[[j]]<-Matrix(W[[j]]%*%X,sparse=TRUE)
    }

    XtWtWX<-list()
    for(j in 1:q){
      XtWtWX[[j]]<-Matrix(t(W[[j]]%*%X)%*%W[[j]]%*%X,sparse=TRUE)
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    sigma_sq_beta<-array(NA,dim=c(q,num_chains,num_samples))
    sigma_sq_lambda<-matrix(NA,num_chains,num_samples)
    sigma_sq<-matrix(NA,num_chains,num_samples)
    h_sq<-array(NA,dim=c(q,num_chains,num_samples))
    alpha_mat<-array(NA,dim=c(p,num_chains,num_samples))
    u_mat<-array(NA,dim=c(m,num_chains,num_samples))
    beta_mat<-array(0,dim=c(n,num_chains,q))

    sigma_sq_beta[,1:num_chains,1]<-1/(2+q)*var(y_scaled)
    sigma_sq_lambda[1:num_chains,1]<-1/(2+q)*var(y_scaled)
    sigma_sq[1:num_chains,1]<-1/(2+q)*var(y_scaled)
    for(cc in 1:num_chains){
      numers<-rep(NA,q)
      denom<-rep(NA,q)
      for(j in 1:q){
        numers[j]<-sigma_sq_beta[j,cc,1]*sum(diag(WGWt[[j]]))
      }
      denom<-sum(numers)+n*sigma_sq[cc,1]
      h_sq[,cc,1]<-numers/denom
      alpha_mat[,cc,1]<-colSums(t(ZtZinvZt)*y_scaled)
      u_mat[,cc,1]<-colSums(t(chol2inv(chol(RtR))%*%t(R))*(y_scaled-colSums(t(Z)*alpha_mat[,cc,1])))
    }
    for(cc in 1:num_chains){
      for(j in 1:q){
        Cov_temp<-chol2inv(chol((1/sigma_sq[cc,1])*XtWtWX[[j]]+(1/sigma_sq_beta[j,cc,1])*In))
        mat_sum_not_j<-rep(0,n)
        if(q>1){
            for(k in c(1:q)[-j]){
                mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
            }
        }
        mu_temp<-(1/sigma_sq[cc,1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,1]-R%*%u_mat[,cc,1]-mat_sum_not_j)
        beta_mat[,cc,j]<-as.vector(mu_temp)
      }
    }

    for(i in 2:num_samples){
      for(cc in 1:num_chains){
        for(j in 1:q){
            tau_beta<-rgamma(1,shape=n/2+a_g,rate=1/2*sum(beta_mat[,cc,j]^2)+b_g)
            sigma_sq_beta[j,cc,i]<-1/tau_beta

            Cov_temp<-chol2inv(chol((1/sigma_sq[cc,i-1])*XtWtWX[[j]]+(1/sigma_sq_beta[j,cc,i])*In))
            mat_sum_not_j<-rep(0,n)
            if(q>1){
                for(k in c(1:q)[-j]){
                    mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
                }
            }
            mu_temp<-(1/sigma_sq[cc,i-1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,i-1]-R%*%u_mat[,cc,i-1]-mat_sum_not_j)
            beta_mat[,cc,j]<-as.vector(mu_temp)+as.vector(t(chol(Cov_temp))%*%rnorm(n))
        }

        gntc_eff<-rep(0,n)
        for(j in 1:q){
            gntc_eff<-gntc_eff+WX[[j]]%*%beta_mat[,cc,j]
        }

        tau_lambda<-rgamma(1,shape=m/2+a_r,rate=1/2*sum(u_mat[,cc,i-1]^2)+b_r)
        sigma_sq_lambda[cc,i]<-1/tau_lambda

        tau<-rgamma(1,shape=n/2+a_e,rate=1/2*sum((y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-colSums(t(R)*u_mat[,cc,i-1])-gntc_eff)^2)+b_e)
        sigma_sq[cc,i]<-1/tau

        u_cov<-chol2inv(chol(1/sigma_sq[cc,i]*RtR+diag(1/sigma_sq_lambda[cc,i],m)))
        u_mean<-as.vector(1/sigma_sq[cc,i]*u_cov%*%t(R)%*%(y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-gntc_eff))
        u_mat[,cc,i]<-mvrnorm(1,u_mean,u_cov)

        alpha_cov<-sigma_sq[cc,i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y_scaled-colSums(t(R)*u_mat[,cc,i])-gntc_eff))
        alpha_mat[,cc,i]<-mvrnorm(1,alpha_mean,alpha_cov)

        numers<-rep(NA,q)
        denom<-rep(NA,q)
        for(j in 1:q){
          numers[j]<-sigma_sq_beta[j,cc,i]*sum(diag(WGWt[[j]]))
        }
        denom<-sum(numers)+n*sigma_sq[cc,i]+sum(diag(R%*%diag(sigma_sq_lambda[cc,i],m)%*%t(R)))
        h_sq[,cc,i]<-numers/denom

      }
      if((num_samples-i)%%100==0){
        cat(i,"/", num_samples, "samples.\n")
      }
    }
    alpha_mat_rescale<-alpha_mat*1/scl_fctr
    u_mat_rescale<-u_mat*1/scl_fctr
    sigma_sq_rescale<-sigma_sq*1/scl_fctr^2
    sigma_sq_lambda_rescale<-sigma_sq_lambda*1/scl_fctr^2
    sigma_sq_beta_rescale<-sigma_sq_beta*1/scl_fctr^2

    rownames(alpha_mat_rescale)<-cov_names1
    colnames(alpha_mat_rescale)<-paste("MCMC_chain_",1:num_chains)
    rownames(u_mat_rescale)<-cov_names2
    colnames(u_mat_rescale)<-paste("MCMC_chain_",1:num_chains)


    rownames(sigma_sq_lambda_rescale)<-paste("MCMC_chain_",1:num_chains)

    rownames(sigma_sq_beta_rescale)<-cov_names3
    colnames(sigma_sq_beta_rescale)<-paste("MCMC_chain_",1:num_chains)
    return(list(alpha=alpha_mat_rescale,u=u_mat_rescale,sigma_sq=sigma_sq_rescale,sigma_sq_lambda=sigma_sq_lambda_rescale,sigma_sq_beta=sigma_sq_beta_rescale))
}


MCMC_heritability_mixed<-function(fixed_form,random_form,dat,G,num_chains=3,num_samples=10000,hyperpars=c(1,1)){

    # 12-18-2019: I set default prior on precision parameter to be
    #             Gamma(1,1). The implied prior on heritability is
    #             Uniform([0,1]).

    a<-hyperpars[1]
    b<-hyperpars[2]

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    scl_fctr<-1/sd(y)
    y_scaled<-y*scl_fctr
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]

    

    ZtZinv<-solve(t(Z)%*%Z)
    ZtZinvZt<-ZtZinv%*%t(Z)

    q<-length(W)
    p<-ncol(Z)
    n<-length(y_scaled)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    svdG<-svd(G)

    X<-Matrix(svdG$u%*%diag(sqrt(svdG$d)),sparse=TRUE)

    XtX<-Matrix(t(X)%*%X,sparse=TRUE)
    In<-Matrix(diag(1,n),sparse=TRUE)

    WX<-list()
    for(j in 1:q){
      WX[[j]]<-Matrix(W[[j]]%*%X,sparse=TRUE)
    }

    XtWtWX<-list()
    for(j in 1:q){
      XtWtWX[[j]]<-Matrix(t(W[[j]]%*%X)%*%W[[j]]%*%X,sparse=TRUE)
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    sigma_sq_lambda<-array(NA,dim=c(length(W),num_chains,num_samples))
    sigma_sq<-matrix(NA,num_chains,num_samples)
    h_sq<-array(NA,dim=c(length(W),num_chains,num_samples))
    alpha_mat<-array(NA,dim=c(p,num_chains,num_samples))
    beta_mat<-array(0,dim=c(n,num_chains,q))

    sigma_sq_lambda[,1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    sigma_sq[1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    for(cc in 1:num_chains){
      numers<-rep(NA,q)
      denom<-rep(NA,q)
      for(j in 1:q){
        numers[j]<-sigma_sq_lambda[j,cc,1]*sum(diag(WGWt[[j]]))
      }
      denom<-sum(numers)+n*sigma_sq[cc,1]
      h_sq[,cc,1]<-numers/denom
      alpha_mat[,cc,1]<-colSums(t(ZtZinvZt)*y_scaled)
    }
    for(cc in 1:num_chains){
      for(j in 1:q){
        Cov_temp<-solve((1/sigma_sq[cc,1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,1])*In)
        mat_sum_not_j<-rep(0,n)
        if(q>1){
            for(k in c(1:q)[-j]){
                mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
            }
        }
        mu_temp<-(1/sigma_sq[cc,1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,1]-mat_sum_not_j)
        beta_mat[,cc,j]<-as.vector(mu_temp)
      }
    }

    for(i in 2:num_samples){
      for(cc in 1:num_chains){
        for(j in 1:q){
            tau_lambda<-rgamma(1,shape=n/2+a,rate=1/2*sum(beta_mat[,cc,j]^2)+b)
            sigma_sq_lambda[j,cc,i]<-1/tau_lambda

            Cov_temp<-solve((1/sigma_sq[cc,i-1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,i])*In)
            mat_sum_not_j<-rep(0,n)
            if(q>1){
                for(k in c(1:q)[-j]){
                    mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
                }
            }
            mu_temp<-(1/sigma_sq[cc,i-1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,i-1]-mat_sum_not_j)
            beta_mat[,cc,j]<-as.vector(mu_temp)+as.vector(t(chol(Cov_temp))%*%rnorm(n))
        }

        rndm_eff<-rep(0,n)
        for(j in 1:q){
            rndm_eff<-rndm_eff+WX[[j]]%*%beta_mat[,cc,j]
        }

        tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-rndm_eff)^2)+b)
        sigma_sq[cc,i]<-1/tau

        alpha_cov<-sigma_sq[cc,i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y_scaled-rndm_eff))

        alpha_mat[,cc,i]<-mvrnorm(1,alpha_mean,alpha_cov)

        numers<-rep(NA,q)
        denom<-rep(NA,q)
        for(j in 1:q){
          numers[j]<-sigma_sq_lambda[j,cc,i]*sum(diag(WGWt[[j]]))
        }
        denom<-sum(numers)+n*sigma_sq[cc,i]
        h_sq[,cc,i]<-numers/denom

      }
      if((num_samples-i)%%1000==0){cat(num_samples-i,' samples remaining.\n')}
    }
    alpha_mat_rescale<-alpha_mat*1/scl_fctr
    sigma_sq_rescale<-sigma_sq*1/scl_fctr^2
    sigma_sq_lambda_rescale<-sigma_sq_lambda*1/scl_fctr^2

    rownames(alpha_mat_rescale)<-cov_names1
    colnames(alpha_mat_rescale)<-paste("MCMC_chain_",1:num_chains)
    rownames(sigma_sq_lambda_rescale)<-cov_names2
    colnames(sigma_sq_lambda_rescale)<-paste("MCMC_chain_",1:num_chains)
    return(list(alpha=alpha_mat_rescale,sigma_sq=sigma_sq_rescale,sigma_sq_lambda=sigma_sq_lambda_rescale))
}



MCMC_heritability_mixed_v2<-function(fixed_form,random_form,dat,G,num_chains=3,num_samples=10000,hyperpars=c(1,1)){

    # 12-18-2019: I set default prior on precision parameter to be
    #             Gamma(1,1). The implied prior on heritability is
    #             Uniform([0,1]).

    # 05-09-2020: Re-derived full conditionals -- trying to see if 
    #             there were any mistakes in the derivations.

    a<-hyperpars[1]
    b<-hyperpars[2]

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    scl_fctr<-1/sd(y)
    y_scaled<-y*scl_fctr
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]

    

    ZtZinv<-solve(t(Z)%*%Z)
    ZtZinvZt<-ZtZinv%*%t(Z)

    q<-length(W)
    p<-ncol(Z)
    n<-length(y_scaled)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    svdG<-svd(G)

    In<-Matrix(diag(1,n),sparse=TRUE)

    WGWt<-list()
    WGWtinv<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
      WGWtinv[[j]]<-solve(WGWt[[j]])
    }

    sigma_sq_lambda<-array(NA,dim=c(length(W),num_chains,num_samples))
    sigma_sq<-matrix(NA,num_chains,num_samples)
    h_sq<-array(NA,dim=c(length(W),num_chains,num_samples))
    alpha_mat<-array(NA,dim=c(p,num_chains,num_samples))
    lambda_mat<-array(0,dim=c(n,num_chains,q))

    sigma_sq_lambda[,1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    sigma_sq[1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    for(cc in 1:num_chains){
      numers<-rep(NA,q)
      denom<-rep(NA,q)
      for(j in 1:q){
        numers[j]<-sigma_sq_lambda[j,cc,1]*sum(diag(WGWt[[j]]))
      }
      denom<-sum(numers)+n*sigma_sq[cc,1]
      h_sq[,cc,1]<-numers/denom
      alpha_mat[,cc,1]<-colSums(t(ZtZinvZt)*y_scaled)
    }
    for(cc in 1:num_chains){
      for(j in 1:q){
        Cov_temp<-solve((1/sigma_sq[cc,1])*In+(1/sigma_sq_lambda[j,cc,1])*WGWtinv[[j]])
        mat_sum_not_j<-rep(0,n)
        if(q>1){
            for(k in c(1:q)[-j]){
                mat_sum_not_j<-mat_sum_not_j+lambda_mat[,cc,k]
            }
        }
        mu_temp<-(1/sigma_sq[cc,1])*Cov_temp%*%(y_scaled-Z%*%alpha_mat[,cc,1]-mat_sum_not_j)
        lambda_mat[,cc,j]<-as.vector(mu_temp)
      }
    }

    for(i in 2:num_samples){
      for(cc in 1:num_chains){
        for(j in 1:q){
            tau_lambda<-rgamma(1,shape=n/2+a,rate=1/2*as.vector(t(lambda_mat[,cc,j])%*%WGWtinv[[j]]%*%lambda_mat[,cc,j])+b)
            sigma_sq_lambda[j,cc,i]<-1/tau_lambda

            Cov_temp<-solve((1/sigma_sq[cc,i-1])*In+(1/sigma_sq_lambda[j,cc,i])*WGWtinv[[j]])
            mat_sum_not_j<-rep(0,n)
            if(q>1){
                for(k in c(1:q)[-j]){
                    mat_sum_not_j<-mat_sum_not_j+lambda_mat[,cc,k]
                }
            }
            mu_temp<-(1/sigma_sq[cc,i-1])*Cov_temp%*%(y_scaled-Z%*%alpha_mat[,cc,i-1]-mat_sum_not_j)
            lambda_mat[,cc,j]<-as.vector(mu_temp)+as.vector(t(chol(Cov_temp))%*%rnorm(n))
        }

        rndm_eff<-rep(0,n)
        for(j in 1:q){
            rndm_eff<-rndm_eff+lambda_mat[,cc,j]
        }

        tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-rndm_eff)^2)+b)
        sigma_sq[cc,i]<-1/tau

        alpha_cov<-sigma_sq[cc,i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y_scaled-rndm_eff))

        alpha_mat[,cc,i]<-mvrnorm(1,alpha_mean,alpha_cov)

        numers<-rep(NA,q)
        denom<-rep(NA,q)
        for(j in 1:q){
          numers[j]<-sigma_sq_lambda[j,cc,i]*sum(diag(WGWt[[j]]))
        }
        denom<-sum(numers)+n*sigma_sq[cc,i]
        h_sq[,cc,i]<-numers/denom

      }
      if((num_samples-i)%%1000==0){cat(num_samples-i,' samples remaining.\n')}
    }
    alpha_mat_rescale<-alpha_mat*1/scl_fctr
    sigma_sq_rescale<-sigma_sq*1/scl_fctr^2
    sigma_sq_lambda_rescale<-sigma_sq_lambda*1/scl_fctr^2

    rownames(alpha_mat_rescale)<-cov_names1
    colnames(alpha_mat_rescale)<-paste("MCMC_chain_",1:num_chains)
    rownames(sigma_sq_lambda_rescale)<-cov_names2
    colnames(sigma_sq_lambda_rescale)<-paste("MCMC_chain_",1:num_chains)
    return(list(alpha=alpha_mat_rescale,sigma_sq=sigma_sq_rescale,sigma_sq_lambda=sigma_sq_lambda_rescale))
}





test_MCMC_heritability_mixed<-function(fixed_form,random_form,dat,G,num_chains=5,num_samples=10000){

    # 10-16-2019: To make prior on heritability to be uniform distribution on
    #             [0,1], I change a from -1 to 1 and b from 0 to 0.001

    # 10-17-2019: Changing a from -1 to 1 on 10-16-2019 gave very poor
    #             estimates of heritability and highly autocorrelated
    #             samples in MCMC chain, so I reverted to a = -1.
    #             The estimates seemed to improve.

    # 10-22-2019: I do not scale in this testing function so I can compare trace plots to "true" 
    #             parameters specified in simulated data.


    dev.new(width=15,height=10)

    a<- -1
    b<-0.001

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    scl_fctr<-1
    y_scaled<-y*scl_fctr
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]

    

    ZtZinv<-solve(t(Z)%*%Z)
    ZtZinvZt<-ZtZinv%*%t(Z)

    q<-length(W)
    p<-ncol(Z)
    n<-length(y_scaled)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    svdG<-svd(G)

    X<-Matrix(svdG$u%*%diag(sqrt(svdG$d)),sparse=TRUE)

    XtX<-Matrix(t(X)%*%X,sparse=TRUE)
    In<-Matrix(diag(1,n),sparse=TRUE)

    WX<-list()
    for(j in 1:q){
      WX[[j]]<-Matrix(W[[j]]%*%X,sparse=TRUE)
    }

    XtWtWX<-list()
    for(j in 1:q){
      XtWtWX[[j]]<-Matrix(t(W[[j]]%*%X)%*%W[[j]]%*%X,sparse=TRUE)
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    sigma_sq_lambda<-array(NA,dim=c(length(W),num_chains,num_samples))
    sigma_sq<-matrix(NA,num_chains,num_samples)
    h_sq<-array(NA,dim=c(length(W),num_chains,num_samples))
    alpha_mat<-array(NA,dim=c(p,num_chains,num_samples))
    beta_mat<-array(0,dim=c(n,num_chains,q))

    sigma_sq_lambda[,1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    sigma_sq[1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    for(cc in 1:num_chains){
      numers<-rep(NA,q)
      denom<-rep(NA,q)
      for(j in 1:q){
        numers[j]<-sigma_sq_lambda[j,cc,1]*sum(diag(WGWt[[j]]))
      }
      denom<-sum(numers)+n*sigma_sq[cc,1]
      h_sq[,cc,1]<-numers/denom
      alpha_mat[,cc,1]<-colSums(t(ZtZinvZt)*y_scaled)
    }
    for(cc in 1:num_chains){
      for(j in 1:q){
        Cov_temp<-solve((1/sigma_sq[cc,1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,1])*In)
        mat_sum_not_j<-rep(0,n)
        if(q>1){
            for(k in c(1:q)[-j]){
                mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
            }
        }
        mu_temp<-(1/sigma_sq[cc,1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,1]-mat_sum_not_j)
        beta_mat[,cc,j]<-as.vector(mu_temp)
      }
    }

    for(i in 2:num_samples){
      for(cc in 1:num_chains){
        for(j in 1:q){
            tau_lambda<-rgamma(1,shape=n/2+a,rate=1/2*sum(beta_mat[,cc,j]^2)+b)
            sigma_sq_lambda[j,cc,i]<-1/tau_lambda

            Cov_temp<-solve((1/sigma_sq[cc,i-1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,i])*In)
            mat_sum_not_j<-rep(0,n)
            if(q>1){
                for(k in c(1:q)[-j]){
                    mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
                }
            }
            mu_temp<-(1/sigma_sq[cc,i-1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,i-1]-mat_sum_not_j)
            beta_mat[,cc,j]<-as.vector(mu_temp)+as.vector(t(chol(Cov_temp))%*%rnorm(n))
        }

        rndm_eff<-rep(0,n)
        for(j in 1:q){
            rndm_eff<-rndm_eff+WX[[j]]%*%beta_mat[,cc,j]
        }

        tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-rndm_eff)^2)+b)
        sigma_sq[cc,i]<-1/tau

        alpha_cov<-sigma_sq[i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y_scaled-rndm_eff))

        alpha_mat[,cc,i]<-mvrnorm(1,alpha_mean,alpha_cov)

        numers<-rep(NA,q)
        denom<-rep(NA,q)
        for(j in 1:q){
          numers[j]<-sigma_sq_lambda[j,cc,i]*sum(diag(WGWt[[j]]))
        }
        denom<-sum(numers)+n*sigma_sq[cc,i]
        h_sq[,cc,i]<-numers/denom

      }
      if((num_samples-i)%%1000==0){cat(num_samples-i,' samples remaining.\n')}
      if((i%%100)==0){
        par(mfrow=c(3,2),fg="white",bg="black",col.axis="white",col.lab="white")
        par2plot<-(sigma_sq_lambda[1,,1:i]+0.25*sigma_sq_lambda[2,,1:i])/(sigma_sq_lambda[1,,1:i]+0.25*sigma_sq_lambda[2,,1:i]+sigma_sq[,1:i])
        plot.new()
        plot.window(xlim=c(0,length(1:i)),ylim=range(par2plot))
        title(main="heritability")
        for(cc in 1:num_chains){
          lines(par2plot[cc,],col=rainbow(num_chains,start=0,end=0.7)[cc],lty=3)
        }
        axis(1)
        axis(2)
        box()

        par2plot<-sigma_sq_lambda[1,,1:i]
        plot.new()
        plot.window(xlim=c(0,length(1:i)),ylim=range(par2plot))
        title(main="sigma_sq_intercept")
        for(cc in 1:num_chains){
          lines(par2plot[cc,],col=rainbow(num_chains,start=0,end=0.7)[cc],lty=3)
        }
        axis(1)
        axis(2)
        box()

        par2plot<-sigma_sq_lambda[2,,1:i]
        plot.new()
        plot.window(xlim=c(0,length(1:i)),ylim=range(par2plot))
        title(main="sigma_sq_z")
        for(cc in 1:num_chains){
          lines(par2plot[cc,],col=rainbow(num_chains,start=0,end=0.7)[cc],lty=3)
        }
        axis(1)
        axis(2)
        box()

        par2plot<-sigma_sq[,1:i]
        plot.new()
        plot.window(xlim=c(0,length(1:i)),ylim=range(par2plot))
        title(main="sigma_sq_error")
        for(cc in 1:num_chains){
          lines(par2plot[cc,],col=rainbow(num_chains,start=0,end=0.7)[cc],lty=3)
        }
        axis(1)
        axis(2)
        box()

        par2plot<-sigma_sq_lambda[1,,1:i]+0.25*sigma_sq_lambda[2,,1:i]+sigma_sq[,1:i]
        plot.new()
        plot.window(xlim=c(0,length(1:i)),ylim=range(par2plot))
        title(main="total variance")
        for(cc in 1:num_chains){
          lines(par2plot[cc,],col=rainbow(num_chains,start=0,end=0.7)[cc],lty=3)
        }
        axis(1)
        axis(2)
        box()
      }
    }
    alpha_mat_rescale<-alpha_mat*1/scl_fctr
    sigma_sq_rescale<-sigma_sq*1/scl_fctr^2
    sigma_sq_lambda_rescale<-sigma_sq_lambda*1/scl_fctr^2

    rownames(alpha_mat_rescale)<-cov_names1
    colnames(alpha_mat_rescale)<-paste("MCMC_chain_",1:num_chains)
    rownames(sigma_sq_lambda_rescale)<-cov_names2
    colnames(sigma_sq_lambda_rescale)<-paste("MCMC_chain_",1:num_chains)
    return(list(alpha=alpha_mat_rescale,sigma_sq=sigma_sq_rescale,sigma_sq_lambda=sigma_sq_lambda_rescale))
}




MCMC_heritability_mixed_with_fixed_pars<-function(fixed_form,random_form,dat,G,num_chains=3,num_samples=10000,hyperpars=c(1,1),parList){

    # 05-07-20: In this function, I request a vector object from the 
    #           user that (optionally) specifies fixed, known values
    #           for the variance parameters. If the values are known,
    #           they will not be sampled. If the values are NA, they
    #           will be sampled.

    a<-hyperpars[1]
    b<-hyperpars[2]

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    scl_fctr<-1/sd(y)
    y_scaled<-y*scl_fctr
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]

    

    ZtZinv<-solve(t(Z)%*%Z)
    ZtZinvZt<-ZtZinv%*%t(Z)

    q<-length(W)
    p<-ncol(Z)
    n<-length(y_scaled)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    svdG<-svd(G)

    X<-Matrix(svdG$u%*%diag(sqrt(svdG$d)),sparse=TRUE)

    XtX<-Matrix(t(X)%*%X,sparse=TRUE)
    In<-Matrix(diag(1,n),sparse=TRUE)

    WX<-list()
    for(j in 1:q){
      WX[[j]]<-Matrix(W[[j]]%*%X,sparse=TRUE)
    }

    XtWtWX<-list()
    for(j in 1:q){
      XtWtWX[[j]]<-Matrix(t(W[[j]]%*%X)%*%W[[j]]%*%X,sparse=TRUE)
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    sigma_sq_lambda<-array(NA,dim=c(length(W),num_chains,num_samples))
    sigma_sq<-matrix(NA,num_chains,num_samples)
    h_sq<-array(NA,dim=c(length(W),num_chains,num_samples))
    alpha_mat<-array(NA,dim=c(p,num_chains,num_samples))
    beta_mat<-array(0,dim=c(n,num_chains,q))

    for(j in 1:q){
      if(is.na(parList[j])){
        sigma_sq_lambda[j,1:num_chains,1]<-(1/2*var(fit_lm1$resid))
      }
      if(!is.na(parList[j])){
        sigma_sq_lambda[j,1:num_chains,]<-parList[j]
      }
    }
  
    if(is.na(parList[q+1])){
      sigma_sq[1:num_chains,1]<-(1/2*var(fit_lm1$resid))
    }
    if(!is.na(parList[q+1])){
      sigma_sq[1:num_chains,]<-parList[q+1]
    }
    for(cc in 1:num_chains){
      numers<-rep(NA,q)
      denom<-rep(NA,q)
      for(j in 1:q){
        numers[j]<-sigma_sq_lambda[j,cc,1]*sum(diag(WGWt[[j]]))
      }
      denom<-sum(numers)+n*sigma_sq[cc,1]
      h_sq[,cc,1]<-numers/denom
      alpha_mat[,cc,1]<-colSums(t(ZtZinvZt)*y_scaled)
    }
    for(cc in 1:num_chains){
      for(j in 1:q){
        Cov_temp<-solve((1/sigma_sq[cc,1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,1])*In)
        mat_sum_not_j<-rep(0,n)
        if(q>1){
            for(k in c(1:q)[-j]){
                mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
            }
        }
        mu_temp<-(1/sigma_sq[cc,1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,1]-mat_sum_not_j)
        beta_mat[,cc,j]<-as.vector(mu_temp)
      }
    }

    for(i in 2:num_samples){
      for(cc in 1:num_chains){
        for(j in 1:q){
            if(is.na(parList[j])){
              tau_lambda<-rgamma(1,shape=n/2+a,rate=1/2*sum(beta_mat[,cc,j]^2)+b)
              sigma_sq_lambda[j,cc,i]<-1/tau_lambda
            }

            Cov_temp<-solve((1/sigma_sq[cc,i-1])*XtWtWX[[j]]+(1/sigma_sq_lambda[j,cc,i])*In)
            mat_sum_not_j<-rep(0,n)
            if(q>1){
                for(k in c(1:q)[-j]){
                    mat_sum_not_j<-mat_sum_not_j+WX[[k]]%*%beta_mat[,cc,k]
                }
            }
            mu_temp<-(1/sigma_sq[cc,i-1])*Cov_temp%*%t(WX[[j]])%*%(y_scaled-Z%*%alpha_mat[,cc,i-1]-mat_sum_not_j)
            beta_mat[,cc,j]<-as.vector(mu_temp)+as.vector(t(chol(Cov_temp))%*%rnorm(n))
        }

        rndm_eff<-rep(0,n)
        for(j in 1:q){
            rndm_eff<-rndm_eff+WX[[j]]%*%beta_mat[,cc,j]
        }

        if(is.na(parList[q+1])){
          tau<-rgamma(1,shape=n/2+a,rate=1/2*sum((y_scaled-colSums(t(Z)*alpha_mat[,cc,i-1])-rndm_eff)^2)+b)
          sigma_sq[cc,i]<-1/tau
        }

        alpha_cov<-sigma_sq[i]*ZtZinv
        alpha_mean<-colSums(t(ZtZinvZt)*(y_scaled-rndm_eff))

        alpha_mat[,cc,i]<-mvrnorm(1,alpha_mean,alpha_cov)

        numers<-rep(NA,q)
        denom<-rep(NA,q)
        for(j in 1:q){
          numers[j]<-sigma_sq_lambda[j,cc,i]*sum(diag(WGWt[[j]]))
        }
        denom<-sum(numers)+n*sigma_sq[cc,i]
        h_sq[,cc,i]<-numers/denom

      }
      if((num_samples-i)%%1000==0){cat(num_samples-i,' samples remaining.\n')}
    }

    alpha_mat_rescale<-alpha_mat*1/scl_fctr

    sigma_sq_rescale<-sigma_sq
    if(is.na(parList[q+1])){
      sigma_sq_rescale<-sigma_sq*1/scl_fctr^2
    }
    sigma_sq_lambda_rescale<-sigma_sq_lambda
    for(j in 1:q){
      if(is.na(parList[j])){
        sigma_sq_lambda_rescale[j,,]<-sigma_sq_lambda[j,,]*1/scl_fctr^2
      }
    }
    rownames(alpha_mat_rescale)<-cov_names1
    colnames(alpha_mat_rescale)<-paste("MCMC_chain_",1:num_chains)
    rownames(sigma_sq_lambda_rescale)<-cov_names2
    colnames(sigma_sq_lambda_rescale)<-paste("MCMC_chain_",1:num_chains)
    return(list(alpha=alpha_mat_rescale,sigma_sq=sigma_sq_rescale,sigma_sq_lambda=sigma_sq_lambda_rescale))
}



bootstrap_heritability<-function(phenotype,G,num_bs=1000){

    y<-phenotype
    n<-length(y)
 
    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    In<-diag(1,n)

    sigma_sq_gamma<-rep(NA,num_bs)
    sigma_sq<-rep(NA,num_bs)
    h_sq<-rep(NA,num_bs)
    mu<-rep(NA,num_bs)

    regress_fit<-regress(y~1,~G)
    sigma_sq_gamma_hat<-regress_fit$sigma[1]
    sigma_sq_hat<-regress_fit$sigma[2]
    numer<-sum(diag(sigma_sq_gamma_hat*G))
    denom<-numer+n*sigma_sq_hat
    h_sq_hat<-numer/denom
    mu_hat<-as.vector(regress_fit$beta)

    y_bs<-mvrnorm(num_bs,mu_hat*rep(1,n),sigma_sq_gamma_hat*G+sigma_sq_hat*In)

    for(i in 1:num_bs){
        regress_fit_bs<-regress(y_bs[i,]~1,~G)
        sigma_sq_gamma[i]<-regress_fit_bs$sigma[1]
        sigma_sq[i]<-regress_fit_bs$sigma[2]
        numer<-sum(diag(sigma_sq_gamma[i]*G))
        denom<-numer+n*sigma_sq[i]
        h_sq[i]<-numer/denom
        mu[i]<-regress_fit_bs$beta
        if((num_bs-i)%%100==0){cat(num_bs-i,' samples remaining.\n')}
    }
    return(list(mu_hat=mu_hat,sigma_sq_hat=sigma_sq_hat,sigma_sq_gamma_hat=sigma_sq_gamma_hat,h_sq_hat=h_sq_hat,mu=mu,sigma_sq=sigma_sq,sigma_sq_gamma=sigma_sq_gamma,h_sq=h_sq))
}


heritability_estimate_MCMC<-function(x,FPR=0.05){
  temp<-c(median(x$h_sq),sd(x$h_sq),quantile(x$h_sq,FPR/2),quantile(x$h_sq,1-FPR/2))
  names(temp)<-c('Posterior median','Posterior SD','Lower bound','Upper bound')
  print(temp)
}

heritability_regression_estimate_MCMC<-function(x,FPR=0.05){
  h_est<-c(median(x$h_sq),sd(x$h_sq),quantile(x$h_sq,FPR/2),quantile(x$h_sq,1-FPR/2))
  names(h_est)<-c('Posterior median','Posterior SD','Lower bound','Upper bound')

  reg_est<-cbind(apply(x$alpha,1,median),apply(x$alpha,1,sd),apply(x$alpha,1,function(y){quantile(y,FPR/2)}),apply(x$alpha,1,function(y){quantile(y,1-FPR/2)}))
  colnames(reg_est)<-c('Posterior median','Posterior SD','Lower bound','Upper bound')

  cat('\n\nHeritability estimates\n\n')
  print(h_est)
  cat('\n\nRegression coefficient estimates\n\n')
  print(reg_est)
}

heritability_mixed_estimate_MCMC<-function(x,FPR=0.05){
  h_est<-cbind(rowMeans(x$h_sq),apply(x$h_sq,1,median),apply(x$h_sq,1,sd),apply(x$h_sq,1,function(ii){quantile(ii,FPR/2)}),apply(x$h_sq,1,function(ii){quantile(ii,1-FPR/2)}))
  colnames(h_est)<-c('Posterior mean','Posterior median','Posterior SD','Lower bound','Upper bound')

  sigma_sq_est<-rbind(c(mean(x$sigma_sq),median(x$sigma_sq),sd(x$sigma_sq),quantile(x$sigma_sq,FPR/2),quantile(x$sigma_sq,1-FPR/2)),cbind(rowMeans(x$sigma_sq_lambda),apply(x$sigma_sq_lambda,1,median),apply(x$sigma_sq_lambda,1,sd),apply(x$sigma_sq_lambda,1,function(ii){quantile(ii,FPR/2)}),apply(x$sigma_sq_lambda,1,function(ii){quantile(ii,1-FPR/2)})))
  colnames(sigma_sq_est)<-c('Posterior mean','Posterior median','Posterior SD','Lower bound','Upper bound')
  rownames(sigma_sq_est)[1]<-"Error"

  reg_est<-cbind(rowMeans(x$alpha),apply(x$alpha,1,median),apply(x$alpha,1,sd),apply(x$alpha,1,function(y){quantile(y,FPR/2)}),apply(x$alpha,1,function(y){quantile(y,1-FPR/2)}))
  colnames(reg_est)<-c('Posterior mean','Posterior median','Posterior SD','Lower bound','Upper bound')

  cat('\n\nHeritability estimates\n\n')
  print(h_est)
  cat('\n\nVariance estimates\n\n')
  print(sigma_sq_est)
  cat('\n\nRegression coefficient estimates\n\n')
  print(reg_est)
}


heritability_estimate_bootstrap<-function(x,FPR=0.05){
    temp<-c(x$h_sq_hat,sd(x$h_sq),quantile(x$h_sq,FPR/2),quantile(x$h_sq,1-FPR/2))
    names(temp)<-c('Estimate','Bootstrap SD','Lower bound','Upper bound')
    print(temp)
}


plot_MCMC_results<-function(x){
    dev.new(width=10,height=1/2*10)
    par(mfrow=c(2,3))
    plot(x$sigma_sq,type='l',main='Error variance',ylab='Error variance',xlab='MCMC iteration')
    plot(x$sigma_sq_gamma,type='l',main='Genetic variance',ylab='Genetic variance',xlab='MCMC iteration')
    plot(x$h_sq,type='l',main='Heritability',ylab='Heritability',xlab='MCMC iteration')
    hist(x$sigma_sq,nclass=50,main='Error variance',xlab='Error variance')
    hist(x$sigma_sq_gamma,nclass=50,main='Genetic variance',xlab='Genetic variance')
    hist(x$h_sq,nclass=50,main='Heritability',xlab='Heritability')
}


plot_MCMC_mixed_results<-function(x){
    num_vc<-nrow(x$h_sq)+1
    dev.new(width=10,height=5)
    par(mfrow=c(ceiling(sqrt(num_vc)),ceiling(sqrt(num_vc))))
    plot(x$sigma_sq,type='l',main='Error variance',ylab='Variance',xlab='MCMC iteration')
    for(j in 1:nrow(x$sigma_sq_lambda)){
      plot(x$sigma_sq_lambda[j,],type='l',main=paste('Genetic variance :',rownames(x$sigma_sq_lambda)[j]),ylab='Variance',xlab='MCMC iteration')
    }
    dev.new(width=10,height=5)
    par(mfrow=c(ceiling(sqrt(num_vc-1)),ceiling(sqrt(num_vc-1))))
    for(j in 1:nrow(x$h_sq)){
      plot(x$h_sq[j,],type='l',main=paste('Heritability :',rownames(x$h_sq)[j]),ylab='Heritability',xlab='MCMC iteration',ylim=c(0,1))
    }
    dev.new(width=10,height=5)
    par(mfrow=c(ceiling(sqrt(num_vc)),ceiling(sqrt(num_vc))))
    hist(x$sigma_sq,main='Error variance',xlab='Variance',nclass=50)
    for(j in 1:nrow(x$sigma_sq_lambda)){
      hist(x$sigma_sq_lambda[j,],main=paste('Genetic variance :',rownames(x$sigma_sq_lambda)[j]),xlab='Variance',nclass=50)
    }
    dev.new(width=10,height=5)
    par(mfrow=c(ceiling(sqrt(num_vc-1)),ceiling(sqrt(num_vc-1))))
    for(j in 1:nrow(x$h_sq)){
      hist(x$h_sq[j,],main=paste('Heritability :',rownames(x$h_sq)[j]),xlab='Heritability',nclass=50,xlim=c(0,1))
    }
}


plot_bootstrap_results<-function(x){
    dev.new(width=10,height=1/4*10)
    par(mfrow=c(1,4))
    hist(x$mu,nclass=50,main='Mean pheno.',xlab='Mean pheno.')
    hist(x$sigma_sq,nclass=50,main='Error variance',xlab='Error variance')
    hist(x$sigma_sq_gamma,nclass=50,main='Genetic variance',xlab='Genetic variance')
    hist(x$h_sq,nclass=50,main='Heritability',xlab='Heritability')
}



ML_heritability<-function(fixed_form,random_form,dat,G){

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]


    q<-length(W)
    p<-ncol(Z)
    n<-length(y)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    negLLfunc<-function(parVal,respVar,covMat,M){
      numObs<-dim(M[[1]])[1]
      V_gamma<-matrix(0,numObs,numObs)
      for(j in 1:q){
        V_gamma<-V_gamma+exp(parVal[j])*M[[j]]
      }
      V_gamma<-V_gamma+exp(parVal[q+1])*diag(1,numObs)
      V_gamma_inv<-solve(V_gamma)
      alpha_gamma<-solve(t(covMat)%*%V_gamma_inv%*%covMat)%*%t(covMat)%*%V_gamma_inv%*%matrix(respVar,numObs,1)

      residVal<- respVar-covMat%*%alpha_gamma
      negLLval<-as.vector(determinant(V_gamma,logarithm=TRUE)$modulus+t(residVal)%*%V_gamma_inv%*%residVal)
      #print(negLLval)
      #print(parVal)
      return(negLLval)
    }

    temp<-list()
    for(i in 1:5){
      temp[[i]]<-optim(rnorm(q+1), negLLfunc, respVar=y, covMat=Z, M=WGWt)
    }
    optimResult<-temp[[which.min(unlist(lapply(temp,function(x){x$value})))]]

    varComp<-optimResult$par
    convCode<-optimResult$convergence
    optVal<-optimResult$value

    V_gamma<-matrix(0,n,n)
    for(j in 1:q){
      V_gamma<-V_gamma+exp(varComp[j])*WGWt[[j]]
    }
    V_gamma<-V_gamma+exp(varComp[q+1])*diag(1,n)
    V_gamma_inv<-solve(V_gamma)
    alpha_gamma<-as.vector(solve(t(Z)%*%V_gamma_inv%*%Z)%*%t(Z)%*%V_gamma_inv%*%matrix(y,n,1))
    sigma_sq_lambda<-exp(varComp[-(q+1)])
    sigma_sq<-exp(varComp[q+1])

    names(alpha_gamma)<-cov_names1
    names(sigma_sq_lambda)<-cov_names2

    return(list(alpha=alpha_gamma,sigma_sq=sigma_sq,sigma_sq_lambda=sigma_sq_lambda,convCode=convCode,optVal=optVal))
}


ML_heritability_v2<-function(fixed_form,random_form,dat,G){

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)

    y<-as.vector(model.response(model.frame(fit_lm1)))
    Z<-model.matrix(fit_lm1)
    W_temp<-model.matrix(fit_lm2)
    W<-list()
    for(j in 1:ncol(W_temp)){
      W[[j]]<-Matrix(diag(W_temp[,j]),sparse=TRUE)
    }
    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]


    q<-length(W)
    p<-ncol(Z)
    n<-length(y)

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
        G<-2*G
    }

    WGWt<-list()
    for(j in 1:q){
      WGWt[[j]]<-Matrix(W[[j]]%*%G%*%t(W[[j]]),sparse=TRUE)
    }

    if(q>2){
      cat("Cannot handle more than 2 random effects.")
    }

    if(q==1){
      hsq_seq<-seq(0,1,0.005)
      ll_seq<-rep(NA,length(hsq_seq))
      for(k in 1:length(hsq_seq)){
        V<-hsq_seq[k]*WGWt[[1]]+(1-hsq_seq[k])*diag(1,n)
        Vinv<-solve(V)
        beta_hat<-solve(t(Z)%*%Vinv%*%Z)%*%t(Z)%*%Vinv%*%y
        tau_sq_hat<-as.vector(t(y-Z%*%beta_hat)%*%Vinv%*%(y-Z%*%beta_hat)/n)

        ll_seq[k]<- -(1/2)*determinant(V*tau_sq_hat)$modulus-(1/2)*t(y-Z%*%beta_hat)%*%(Vinv/tau_sq_hat)%*%(y-Z%*%beta_hat)
      }

      plot(ll_seq~hsq_seq,cex=0.5,xlab="heritability_overall",ylab="log-likelihood")

      hsq_hat<-hsq_seq[which.max(ll_seq)]
      V<-hsq_hat*WGWt[[1]]+(1-hsq_hat)*diag(1,n)
      Vinv<-solve(V)

      beta_hat<-as.vector(solve(t(Z)%*%Vinv%*%Z)%*%t(Z)%*%Vinv%*%y)
      tau_sq_hat<-as.vector(t(y-Z%*%beta_hat)%*%Vinv%*%(y-Z%*%beta_hat)/n)

      names(beta_hat)<-cov_names1

      return(list(regression_coefficients=beta_hat,heritability_overall=hsq_hat,variance_total=tau_sq_hat,hsq_vals=hsq_seq,loglike_vals=ll_seq))
    }

    if(q==2){
      hsq_a_seq<-seq(0,1,0.02)
      hsq_b_seq<-seq(0,1,0.02)
      ll_seq<-matrix(NA,length(hsq_a_seq),length(hsq_b_seq))
      for(j in 1:length(hsq_a_seq)){
        for(k in 1:length(hsq_b_seq)){
          if((hsq_a_seq[j]+hsq_b_seq[k])<=1){
            V<-hsq_a_seq[j]*WGWt[[1]]+hsq_b_seq[k]*WGWt[[2]]+(1-hsq_a_seq[j]-hsq_b_seq[k])*diag(1,n)
            Vinv<-solve(V)
            beta_hat<-solve(t(Z)%*%Vinv%*%Z)%*%t(Z)%*%Vinv%*%y
            tau_sq_hat<-as.vector(t(y-Z%*%beta_hat)%*%Vinv%*%(y-Z%*%beta_hat)/n)

            ll_seq[j,k]<- as.vector(-(1/2)*determinant(V*tau_sq_hat)$modulus-(1/2)*t(y-Z%*%beta_hat)%*%(Vinv/tau_sq_hat)%*%(y-Z%*%beta_hat))
          }
          if((hsq_a_seq[j]+hsq_b_seq[k])>1){
            ll_seq[j,k]<- -Inf
          }
        }
      }
      rownames(ll_seq)<-hsq_a_seq
      colnames(ll_seq)<-hsq_b_seq

      ll_seq[ll_seq==-Inf]<-min(ll_seq[ll_seq>-Inf])

      heatmap(t(ll_seq),Rowv=NA,Colv=NA,scale="none",xlab="heritability_overall",ylab="heritability_regression")

      hsq_a_hat<-as.numeric(names(which.max(apply(ll_seq,1,max))))   # hsq_a_hat
      hsq_b_hat<-as.numeric(names(which.max(apply(ll_seq,2,max))))   # hsq_b_hat

      V<-hsq_a_hat*WGWt[[1]]+hsq_b_hat*WGWt[[2]]+(1-hsq_a_hat-hsq_b_hat)*diag(1,n)
      Vinv<-solve(V)

      beta_hat<-as.vector(solve(t(Z)%*%Vinv%*%Z)%*%t(Z)%*%Vinv%*%y)
      tau_sq_hat<-as.vector(t(y-Z%*%beta_hat)%*%Vinv%*%(y-Z%*%beta_hat)/n)

      names(beta_hat)<-cov_names1

      return(list(regression_coefficients=beta_hat,heritability_overall=hsq_a_hat,heritability_regression=hsq_b_hat,variance_total=tau_sq_hat,hsq_a_vals=hsq_a_seq,hsq_b_vals=hsq_b_seq,loglike_vals=ll_seq))
    }
}



ML_heritability_bootstrap_test<-function(fixed_form,random_form_full,random_form_reduced=NULL,dat,G,nBootstraps=100){

    MLest<-ML_heritability(fixed_form,random_form_full,dat,G)

    form1<-fixed_form
    form2<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form_full)[2]))
    form3<-as.formula(paste(as.character(fixed_form)[2],"~",as.character(random_form_reduced)[2]))

    fit_lm1<-lm(form1,data=dat)
    fit_lm2<-lm(form2,data=dat)
    fit_lm3<-lm(form3,data=dat)

    Z<-model.matrix(fit_lm1)
    W_temp_reduced<-model.matrix(fit_lm3)

    if(ncol(W_temp_reduced)==0){
      W_reduced<-NULL
    }
    if(ncol(W_temp_reduced)>0){
      W_reduced<-list()
      for(j in 1:ncol(W_temp_reduced)){
        W_reduced[[j]]<-Matrix(diag(W_temp_reduced[,j]),sparse=TRUE)
      }
    }

    cov_names1<-dimnames(model.matrix(fit_lm1))[[2]]
    cov_names2<-dimnames(model.matrix(fit_lm2))[[2]]
    cov_names3<-dimnames(model.matrix(fit_lm3))[[2]]


    q_reduced<-length(W_reduced)
    p<-ncol(Z)
    n<-dim(Z)[1]

    if(!is.null(fit_lm1$na.action)){
      G<-G[-fit_lm1$na.action,-fit_lm1$na.action]
    }

    if(abs(mean(diag(G))-0.5)<abs(mean(diag(G))-1)){
      G<-2*G
    }

    if(q_reduced>0){
      WGWt_reduced<-list()
      for(j in 1:q_reduced){
        WGWt_reduced[[j]]<-Matrix(W_reduced[[j]]%*%G%*%t(W_reduced[[j]]),sparse=TRUE)
      }
    }

    V_mat<-matrix(0,n,n)
    if(q_reduced>0){
      for(j in 1:q_reduced){
        parInd<-which(names(MLest$sigma_sq_lambda)==cov_names3[j])
        V_mat<-V_mat+MLest$sigma_sq_lambda[parInd]*WGWt_reduced[[j]]
      }
    }
    V_mat<-V_mat+MLest$sigma_sq*diag(1,n)
    V_mat_inv<-solve(V_mat)
    alphaEst<-MLest$alpha

    MLest_bs<-list()
    for(i in 1:nBootstraps){
      dat$y_bs<-mvrnorm(1,Z%*%alphaEst,V_mat)
      temp<-as.character(fixed_form)
      fixed_form_bs<-as.formula(paste("y_bs","~",temp[3]))
      MLest_bs[[i]]<-ML_heritability(fixed_form_bs,random_form_full,dat,G)
      cat(i,"bootstraps finished.\n")
    }

    components2test<-cov_names2[!(cov_names2%in%cov_names3)]

    obsEst<-rep(NA,length(components2test))   
    bsEst<-matrix(NA,nBootstraps,length(obsEst))
    par(mfrow=c(1,length(obsEst)))
    for(k in 1:length(obsEst)){
      obsEst[k]<-MLest$sigma_sq_lambda[names(MLest$sigma_sq_lambda)==components2test[k]]
      bsEst[,k]<-unlist(lapply(MLest_bs,function(x){x$sigma_sq_lambda[names(x$sigma_sq_lambda)==components2test[k]]}))
      hist(bsEst[,k],nclass=nBootstraps/2,xlim=range(bsEst[,k],obsEst[k]),xlab=paste("Bootstrap distribution,",components2test[k]),main=paste("Bootstrap distribution,",components2test[k]))
      abline(v=obsEst[k],col="firebrick3",lwd=2)
    }
    names(obsEst)<-components2test
    colnames(bsEst)<-names(obsEst)

    if(length(components2test)==1){
      bootstrap_pvalue<-(sum(t(bsEst)>=obsEst)+1)/(nBootstraps+1)    
    }
    if(length(components2test)>1){
      bootstrap_pvalue<-(rowSums(t(bsEst)>=obsEst)+1)/(nBootstraps+1)    
    }

    return(list(pvalue=bootstrap_pvalue,bootstrap_estimates=MLest_bs,ML_estimates=MLest))
}


