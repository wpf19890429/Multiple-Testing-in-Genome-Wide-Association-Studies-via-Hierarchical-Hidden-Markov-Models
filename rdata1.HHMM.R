rdata1.HHMM<-function(NUM, S, pii_eta, pii_theta, A, B, f0, f1)
{

#############################

## Parameters
 # NUM: the number of hypotheses
 # S: the size of the block
 # pii_eta: the initial probabilities of eta
 # pii_theta: the initial probabilities of theta, K*1
 # A: the transition matrix of theta
 # B: the block-wide transition matrix of eta 
 # f0: the pdf of z_i given theta_i=0
 # f1: the pdf of z_i given theta_i=1


## Values
 # theta: the states of the null hypotheses
 # eta: the class of region
 # z: the observations

###############################

## NUM is a multiple of block_size
	eta<-rep(0, NUM)
	theta<-rep(0, NUM)
	z<-rep(0, NUM)


## generating the states

 # initial state
      K<-length(pii_eta)
      eta[1:S]<-rep(sample(1:K, 1, replace=FALSE, prob=pii_eta), S)
      theta[1]<-rbinom(1, 1, prob=pii_theta[eta[1]])

 # other states

      NUM_b<-NUM/S
	for(i in 2:NUM_b)
	{
		block_start<-(i-1)*S+1
            block_end<-i*S
		eta[block_start:block_end]<-rep(sample(1:K, 1, replace=FALSE, prob=B[eta[block_start-1], ]), S)	
	}
	for(i in 2:NUM)
	{
		if(theta[i-1]==0)
			theta[i]<-rbinom(1, 1, A[1, 2, eta[i]])
		else
			theta[i]<-rbinom(1, 1, A[2, 2, eta[i]])
	}



## generating the observations
	for (i in 1:NUM)
	{
		if (theta[i]==0)
			z[i]<-rnorm(1, mean=f0[1], sd=f0[2])
		else
			z[i]<-rnorm(1, mean=f1[1], sd=f1[2])
	}

	data<-list(theta=theta, eta=eta, z=z)
	return(data)

}





