em.HHMM<-function(z, K, S, maxiter=200)
{

#############################

## Usage
 # em.HHMM(z, maxiter=200)


## Arguments

 # z=(z[1], ..., z[m]): the observations
 # maxiter: the maximum number of iterations


## Values

 # niter: number of iterations


#############################

	NUM<-length(z)
# precision tolerance level
	ptol<-1
	niter<-0

### initializing model parameters

	pii_eta.new<-rep(0.5, K)
      pii_theta.new<-rep(0.5, K)
	A.new<-array(0, c(2, 2, K))
	A.new[, ,1]<-matrix(c(0.6, 0.4, 0.4, 0.6), 2, 2, byrow=TRUE)
	A.new[, ,2]<-matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2, byrow=TRUE)
	B.new<-matrix(1/K, nrow=K, ncol=K, byrow=TRUE)
	f0<-c(0, 1)
      f0z<-dnorm(z, mean=0, sd=1)
	f1z.new<-dnorm(z, mean=2, sd=1)


### initializing intermediate variables

	phi_1<-rep(0, K)
	xi<-array(0, c(NUM-1, 2, 2, K, K))
      zeta<-array(0, c(NUM-1, 2, 2, K))
      psi<-array(0, c(NUM-1, 2, K))
	rho_1<-matrix(0, nrow=2, ncol=K)
	gamma<-matrix(0, nrow=NUM, ncol=2)
	nu<-array(0, c(NUM-1, K, K))
	tau<-matrix(0, nrow=NUM-1, ncol=K)

	Loglikelihood.new<--10000
	diff<-100
	ptol<-1e-5


### The E-M Algorithm


	pb<-PB(Methods="em.HHMM", Rep=maxiter)

	while(niter<maxiter)
	{

		niter<-niter+1
	      pb$tick()


		pii_eta.old<-pii_eta.new
      	pii_theta.old<-pii_theta.new
		A.old<-A.new
		B.old<-B.new
		f1z.old<-f1z.new
		Loglikelihood.old<-Loglikelihood.new


## updating the weights and probabilities of hidden states

		bwfw.res<-bwfw.HHMM(z, S, pii_eta.old, pii_theta.old, A.old, B.old, f0, f1z.old)


# the backward-forward variable

		alpha<-bwfw.res$bw
		beta<-bwfw.res$fw


# calculating xi

		b1<-rep(0, NUM-1)
		for(j in 1:(NUM-1))
		{
			if(j%%S==0)
                  {
                  	for(p in 1:2)
                        	for(k in 1:K)
						for(l in 1:K)
                                    {
							xi[j, p, 1, k, l]<-alpha[j, p, k]*f0z[j+1]*beta[j+1, 1, l]*A.old[p, 1, k]*B.old[k, l]
                                          xi[j, p, 2, k, l]<-alpha[j, p, k]*f1z.old[j+1]*beta[j+1, 2, l]*A.old[p, 2, k]*B.old[k, l]
                                          b1[j]<-b1[j]+xi[j, p, 1, k, l]+xi[j, p, 2, k, l]
                                     } 
                  }
			else
                  {
                  	for(p in 1:2)
                        	for(k in 1:K)
                              {
                              	xi[j, p, 1, k, k]<-alpha[j, p, k]*f0z[j+1]*beta[j+1, 1, k]*A.old[p, 1, k]   
						xi[j, p, 2, k, k]<-alpha[j, p, k]*f1z.old[j+1]*beta[j+1, 2, k]*A.old[p, 2, k]   	  
                                    b1[j]<-b1[j]+xi[j, p, 1, k, k]+xi[j, p, 2, k, k]
					}
			}
                  for(p in 1:2)
                  	for(q in 1:2)
					for(k in 1:2)
						for(l in 1:2)
							xi[j, p, q, k, l]<-1/b1[j]*xi[j, p, q, k, l]
		}  


# calculating phi_1

		for(k in 1:K)
			phi_1[k]<-sum(xi[1, 1, 1, k, ])+sum(xi[1, 1, 2, k, ])+sum(xi[1, 2, 1, k, ])+sum(xi[1, 2, 2, k, ])


# calculating zeta

		for(j in 1:(NUM-1))
			for(p in 1:2)
				for(q in 1:2)
					for(l in 1:K)
						zeta[j, p, q, l]<-sum(xi[j, p, q, , l])



# calculating psi

		for(j in 1:(NUM-1))
			for(p in 1:2)
				for(l in 1:K)
					psi[j, p, l]<-sum(zeta[j, p, , l])


# calculating mu

		for(j in 1:(NUM-1))
			for(k in 1:K)
				for(l in 1:K)
					nu[j, k, l]<-xi[j, 1, 1, k, l]+xi[j, 1, 2, k, l]+xi[j, 2, 1, k, l]+xi[j, 2, 2, k, l]	

# calculating tau

		for(j in 1:(NUM-1))
			for(k in 1:K)
				tau[j, k]<-sum(nu[j, k, ])

                          	
# calculating rho_1

		for(p in 1:2)
			for(k in 1:K)
                  	rho_1[p, k]<-sum(xi[1, p, 1, k, ])+sum(xi[1, p, 2, k, ]) 


# calculating gamma 

		gamma[, 1]<-bwfw.res$HLIS
            gamma[, 2]<-1-gamma[, 1]
				
## updating parameters

		pii_eta.new<-phi_1
		pii_theta.new<-rho_1[2, ]/apply(rho_1, 2, sum)
		for(p in 1:2)
			for(q in 1:2)
				for(l in 1:K)
					A.new[p, q, l]<-sum(zeta[, p, q, l])/sum(psi[, p, l])

		index<-which(seq(1, NUM-1, 1)%%S==0)		
		for(k in 1:K)
			for(l in 1:K)
				B.new[k, l]<-sum(nu[index, k, l])/sum(tau[index, k])
                                       
		kern_f1<-density(z, weights=gamma[, 2]/sum(gamma[, 2]))
		f1z.new<-approx(kern_f1$x, kern_f1$y, z, rule=2, ties="ordered")$y


		Loglikelihood.new<-sum(log(bwfw.res$c0))
            if(abs(Loglikelihood.new-Loglikelihood.old)/abs(Loglikelihood.old)<1e-5)
			break


	}

	em.var<-list(pii_eta=pii_eta.new, pii_theta=pii_theta.new, A=A.new, B=B.new, f1z=f1z.new, ni=niter, gamma=gamma)
	return (em.var)

}

