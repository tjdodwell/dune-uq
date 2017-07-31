#ifndef common_stats_h
#define common_stats_h

double inline computeMean(std::vector<double>& samples){
	int N = samples.size();
	double mu = 0.0;
	for (int i = 0; i < N; i++){	mu += samples[i];	}
	mu /= N;
	return mu;
}

double inline computeVar(std::vector<double>& samples){	
	int N = samples.size();
	double mu = 0.0;
	for (int i = 0; i < N; i++){
		mu += samples[i];
	}
	mu /= N;
	double var = 0.0;
	for (int i = 0; i < N; i++){
		var += (samples[i] - mu) * (samples[i] - mu);
	}
	var /= N - 1;
	return var;
}

double inline autocorrelation(std::vector<double>& X, double mu, double s, int k){

	int n = X.size();

	double rho_k = 0.0;

	for (int i = 0; i < n - k; i++){

		rho_k += (X[i] - mu) * (X[i + k] - mu);

	}

	rho_k /= n * s;

	return rho_k;

}

double inline effectiveSampleSize(std::vector<double>& X){

	double tau = 1.0;

	int n = X.size();

	double mu = computeMean(X);

	double s = computeVar(X);

	for (int k = 1; k < n-1; k++){
		tau += 2.0 * autocorrelation(X,mu,s,k);
	}

	return tau;

}

#endif