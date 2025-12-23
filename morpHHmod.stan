data{
	int<lower=0> N; // sample size
	int<lower=0> Npop; // number of populations
	int<lower=0> id[N]; // vector of population ID
	vector[N] y; // response variable, PC score
}

parameters{
	real alpha[Npop]; // population intercepts
	real<lower=0> sig[Npop]; // residual SD
}

model{

	// linear model
	for(i in 1:N){
		y[i] ~ normal(alpha[id[i]], sig[id[i]]);
	}
	// priors
	alpha ~ normal(0,100);	
	sig ~ normal(0,10);	
}

