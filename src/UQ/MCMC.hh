#include "../stats/common.hh"

/* MCMC.hh - Markov Chain Monte Carlos Algorithm */
template <class MODEL, class RandomField, class Samples>
void inline MCMC(MODEL& model, RandomField& z,Samples& samples, int level = 0){

	int rank, nproc;

    MPI_Comm new_comm;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
            std::cout << "----------------------------------" << std::endl;
            std::cout << "-    dune-uq : MCMC		       -" << std::endl;
            std::cout << "----------------------------------" << std::endl;
     }


	int N = config.get<int>("MCMC.initialSamples", 1000);

	double greedyFactor = config.get<double>("MCMC.greedyFactor", 1.1);
	
	double epsilon = config.get<double>("MCMC.samplingError", 1e-6);

	int maxSamples = config.get<int>("MCMC.maxSamples", 10e6);

	int confidence = config.get<double>("MCMC.confidence", 0.9);

	

	Dune::Timer watch;


	// === Burnin Samples

	// Do some initial samples before computing autocorrelation time

	std::vector<double> Q_burnin(N);

	z.generate_random_field(); // Generate random field from prior

	Q_burnin[0] = model.getSample(level,z,false); // z contains likelihood for xi	

	double accept = 0.0;


	for (int i = 1; i < N; i++){

		z.generate_proposal(); // Make proposal 

		double Q_tmp = model.getSample(level,z,true);

		bool accept_flag = z.accept_proposal();

		//z.print_likelihood();

		if(accept_flag){
			Q_burnin[i] = Q_tmp; 
			accept += 1.0;
		}
		else{
			Q_burnin[i] = Q_burnin[i-1];
		}


	}

	double accept_ratio = accept / N;

	std::cout << "Acceptance Ratio = " <<  accept / N << std::endl; 
	std::cout << effectiveSampleSize(Q_burnin) << std::endl;


}

	

