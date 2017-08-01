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

double inline compute_autocorrelation_time(std::vector<double>& X, int m = 30){

double s = computeVar(X); // sample variance

int n = X.size();

std::vector<std::vector<double>> batch(m); // make a vector for each batch

std::vector<double> mu(m); // vector for the means of the batches

int size_of_batch = std::floor(n/m);

for (int j = 0; j < m; j++){ // for each batch

	batch[j].resize(size_of_batch); // resize j^th batch to contain values

	for (int i = 0; i < size_of_batch; i++){
		batch[j][i] = X[j*m + i]; // store values for each batch
	}
	mu[j] = computeMean(batch[j]); // compute mean of jth batch
}

double sm = computeVar(mu); // compute variances

double tau = size_of_batch * sm / s; // estimate autocorrelation time

std::cout << "size_of_batch = " << size_of_batch << std::endl;

std::cout << "sm = " << sm << std::endl;

std::cout << "s = " << s << std::endl;

return tau;
	
}



int inline effectiveSampleSize(std::vector<double>& X, int m = 30){

	int n = X.size();

	double tau = compute_autocorrelation_time(X,m);

	int neff = std::floor(n / tau);

	return neff;

}

#endif