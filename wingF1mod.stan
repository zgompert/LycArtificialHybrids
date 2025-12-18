data{
	int<lower=0> N; // sample size
	int<lower=0> Npop; // number of populations
	int<lower=0> nid1[N]; // vector of population ids for male 
	int<lower=0> nid2[N]; // vector of population ids for female 
	int<lower=0> hid[N]; // vector of boolean hybrid, 1 = yes
	int<lower=0> annaM[N]; // vector of boolean male is anna, 1 = yes
	vector[N] y; // response variable, PC score
}

parameters{
	real alpha[Npop]; // population intercepts
	real betaH; // effect of being a hybrid
	real betaAM; // effect, for hybrids, of anna male
	real sig; // residual SD
}

model{
	real mu[N]; // normal expectations

	// linear model
	for(i in 1:N){
		mu[i] = (alpha[nid1[i]] + alpha[nid2[i]])/2 + hid[i] * betaH +
			hid[i] * annaM[i] * betaAM;
		y[i] ~ normal(mu[i], sig);
	}
	// priors
	alpha ~ normal(0,100);	
	betaH ~ normal(0,20);	
	betaAM ~ normal(0,20);
	sig ~ normal(0,10);	
}

