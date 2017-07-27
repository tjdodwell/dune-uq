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

#endif