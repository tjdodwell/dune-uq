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
            std::cout << "-    dune-uq : Standard MC     -" << std::endl;
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

	std::vector<double> Q_burnin(1000);

	z.generate_random_field(); // Generate random field from prior

	Q_burnin[] = model.getSample(level,z); // z contains likelihood for xi	

	for (int i = 0; i < N; i++){

		z.make_proposal(false); // Make proposal 

		double Q_tmp = model.getSample(level,z);

		if(z.accept_proposal()){
			Q_burnin[i] = Q_tmp; 
		}

	}





	// Do some stats on these intial samples

	double V = computeVar(Q);
	double samplingError = std::sqrt(V/N); // Compute sampling error
	bool flag = false;
	if (samplingError > epsilon){flag = true;}

	int Nold;

	while(flag && N < maxSamples){

		Nold = N;
		N *= greedyFactor; N += 1;

		Q.resize(N);

		for (int i = Nold; i < N; i++){
			watch.reset();
			z.generate_random_field();
			Q[i] = model.getSample(level,z);
			samples.add_Sample(Q[i],z, watch.elapsed());
		}

		V = computeVar(Q);
		samplingError = std::sqrt(V/N);

		if (samplingError < epsilon){flag = false;}

	}

}

	

