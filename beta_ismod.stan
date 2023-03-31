data{
	int L; /* # of SNPs */
	int N; /* # of pops */
	matrix[L, N] p; /* matrix of allele freqs */
}

parameters{
	real<lower=0,upper=50> Nm; /* number of migrants */
	vector<lower=0.0,upper=1.0>[L] pm; /* migrant allele freq. */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment beta for allele freqs. */
			target += beta_lpdf(p[i, j] | 4*Nm*pm[i], 4*Nm*(1-pm[i]));
		}
	}

	for(i in 1:L){
		/* increment beta prior on pm */
		target += beta_lpdf(pm[i] | 0.5, 0.5);
	}

	/* increment cauchy Nm */
	target += cauchy_lpdf(Nm | 0, 10);

}
