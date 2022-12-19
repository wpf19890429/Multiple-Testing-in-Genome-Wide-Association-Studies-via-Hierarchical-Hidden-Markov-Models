bwfw.HHMM<-function(z, S, pii_eta, pii_theta, A, B, f0, f1z)
{

#############################

## Parameters
 # z: the observations
 # S: the size of the block
 # pii_eta: the initial probabilities of eta
 # pii_theta: the initial probabilities of theta, K*2
 # A: the transition matrix of theta
 # B: the block-wide transition matrix of eta 
 # f0: the pdf of z_i given theta_i=0
 # f1z: the probability of z_i given theta_i=1


## Values
 # alpha: rescaled backward variables
 # beta: rescaled forward variables	
 # HLIS: the HLIS statistics

###############################

## Initialize

	NUM<-length(z)
	K<-length(pii_eta)

## Densities

	f0z<-dnorm(z, f0[1], f0[2])


## the backward-forward procedure

# a. the backward variables
# --rescaled 


	alpha<-array(0, dim=c(NUM, 2, K))
# scaling variable c_0
	c0<-rep(0, NUM)

	alpha[1, 1, ]<-pii_eta*(1-pii_theta)*f0z[1]
      alpha[1, 2, ]<-pii_eta*pii_theta*f1z[1]
      c0[1]<-1

	for (i in 1:(NUM-1))
	{ 
		for(k in 1:K)
		{
	            if(i%%S==0)
			{
				alpha[i+1, 1, k]<-f0z[i+1]*(A[1, 1, k]*t(alpha[i, 1, ])%*%B[, k]+A[2, 1, k]*t(alpha[i, 2, ])%*%B[, k])
				alpha[i+1, 2, k]<-f1z[i+1]*(A[1, 2, k]*t(alpha[i, 1, ])%*%B[, k]+A[2, 2, k]*t(alpha[i, 2, ])%*%B[, k])
                  }
                  else
                  {
                        alpha[i+1, 1, k]<-f0z[i+1]*(A[1, 1, k]*alpha[i, 1, k]+A[2, 1, k]*alpha[i, 2, k])
                        alpha[i+1, 2, k]<-f1z[i+1]*(A[1, 2, k]*alpha[i, 1, k]+A[2, 2, k]*alpha[i, 2, k])
                   }
		}
# rescaling alpha
            c0[i+1]<-1/(sum(alpha[i+1, 1, ])+sum(alpha[i+1, 2, ]))            
	      alpha[i+1, 1, ]<-c0[i+1]*alpha[i+1, 1, ]
	      alpha[i+1, 2, ]<-c0[i+1]*alpha[i+1, 2, ]  
	}


# b. the forward variables
# --rescaled

	beta<-array(0, dim=c(NUM, 2, K))
# scaling variable b_0
	b0<-rep(0, NUM)
	beta[NUM, 1, ]<-1/(2*K)
      beta[NUM, 2, ]<-1/(2*K)


	for (i in (NUM-1):1)
	{ 
		for(k in 1:K)
		{
	            if(i%%S==0)
			{
				beta[i, 1, k]<-t(beta[i+1, 1, ])%*%(A[1, 1, ]*B[k, ])*f0z[i+1]+t(beta[i+1, 2, ])%*%(A[1, 2, ]*B[k, ])*f1z[i+1]
				beta[i, 2, k]<-t(beta[i+1, 1, ])%*%(A[2, 1, ]*B[k, ])*f0z[i+1]+t(beta[i+1, 2, ])%*%(A[2, 2, ]*B[k, ])*f1z[i+1]
                  }
                  else
                  {
				beta[i, 1, k]<-beta[i+1, 1, k]*A[1, 1, k]*f0z[i+1]+beta[i+1, 2, k]*A[1, 2, k]*f1z[i+1]
				beta[i, 2, k]<-beta[i+1, 1, k]*A[2, 1, k]*f0z[i+1]+beta[i+1, 2, k]*A[2, 2, k]*f1z[i+1]
                  }
            }
  # rescaling beta
            b0[i]<-1/(sum(beta[i, 1, ])+sum(beta[i, 2, ]))            
	      beta[i, 1, ]<-b0[i]*beta[i, 1, ]
	      beta[i, 2, ]<-b0[i]*beta[i, 2, ]
	} 


	HLIS<-rep(0, NUM)

	for (k in 1:NUM)
	{ 
		q1<-sum(alpha[k, 1, ]*beta[k, 1, ])
		q2<-sum(alpha[k, 2, ]*beta[k, 2, ])
		HLIS[k]<-q1/(q1+q2)
	}

	bwfw.var<-list(bw=alpha, fw=beta, HLIS=HLIS, c0=c0)
	return(bwfw.var)  
}

