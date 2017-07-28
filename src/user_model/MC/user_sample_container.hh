#include "../../stats/common.hh"

template<class RandomField>
class SampleContainer{

public:

	SampleContainer(){

		N = 0;

		stochastic_dim = config.get<int>("RandomField.dim",100);
		

	};

	void inline add_Sample(double newSample, RandomField& z, double t){

		N += 1;

		Q.push_back(newSample);
		time.push_back(t);

		xi.resize(N);

		xi[N-1].resize(stochastic_dim);

		for (int i = 0; i < stochastic_dim; i++){
			xi[N-1][i] = z.getValue(i);
		}

	}

	void inline post(){

		double EQ = computeMean(Q);

		double average_time = computeMean(time);

		double V = computeVar(Q);
		double samplingError = std::sqrt(V/N);

		std::cout << "== Results == "<< std::endl;
		std::cout << "E[Q] = " << EQ << std::endl;
		std::cout << "V[Q] = " << V << std::endl;
		std::cout << "Number of Samples = " << N << std::endl;
		std::cout << "Sampling Error = " << samplingError << std::endl;
		std::cout << "Average Time per Samples = " << average_time << std::endl;

	}

private:

	int N, stochastic_dim;
	std::vector<double> Q, time;
	std::vector<std::vector<double>> xi;
};