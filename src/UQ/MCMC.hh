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

	int burnin_factor = config.get<int>("MCMC.burnin_factor", 10);

	int m = config.get<int>("MCMC.batch_size",30);



	Dune::Timer watch;


	// === Burnin Samples

	// Do some initial samples before computing autocorrelation time for the first time

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

	double tau = compute_autocorrelation_time(Q_burnin,m);

	while (N < std::ceil(burnin_factor * tau)){

		Q_burnin.resize(std::ceil(burnin_factor * tau));

		for (int i = N; i < std::ceil(burnin_factor * tau); i++){

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

		N = std::ceil(burnin_factor * tau);

		tau = compute_autocorrelation_time(Q_burnin,m);

	}

	// Burnin Complete

	int subSamplingRate = 2 * std::ceil(tau);


	std::vector<double> Q(N);

	Q[0] = Q_burnin[N-1];

	double Qsub = 0.0;

	// Do some initial samples
	for (int i = 1; i < N; i++){

		watch.reset();

		Qsub = Q[i-1];

		for (int j = 1; j < subSamplingRate; j++){

			z.generate_proposal(); // Make proposal 

			double Q_tmp = model.getSample(level,z,true);

			bool accept_flag = z.accept_proposal();

			//z.print_likelihood();

			if(accept_flag){
				Qsub = Q_tmp; 
			}
		
		}

		Q[i] = Qsub;

		samples.add_Sample(Q[i],z,watch.elapsed());
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

			Qsub = Q[i-1];

			for (int j = 1; j < subSamplingRate; j++){

				z.generate_proposal(); // Make proposal 

				double Q_tmp = model.getSample(level,z,true);

				bool accept_flag = z.accept_proposal();

				//z.print_likelihood();

				if(accept_flag){
					Qsub = Q_tmp; 
				}
			
			}

			Q[i] = Qsub;

			samples.add_Sample(Q[i],z,watch.elapsed());
		}

		V = computeVar(Q);
		samplingError = std::sqrt(V/N);

		if (samplingError < epsilon){flag = false;}

	}





}

	

