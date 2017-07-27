#include "../stats/common.hh"

/* MC.hh - Simple Monte Carlos Algorithms */
template <class MODEL, class RandomField>
void inline MC(MODEL& model, RandomField& z, int level = 0){

	int rank, nproc;

    MPI_Comm new_comm;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
            std::cout << "----------------------------------" << std::endl;
            std::cout << "-    dune-uq : Standard MC     -" << std::endl;
            std::cout << "----------------------------------" << std::endl;
     }


	int N = config.get<int>("MC.initialSamples", 100);

	double greedyFactor = config.get<double>("MC.greedyFactor", 1.1);
	
	double epsilon = config.get<double>("MC.samplingError", 1e-6);

	int maxSamples = config.get<int>("MC.maxSamples", 10e6);

	int confidence = config.get<double>("MC.confidence", 0.9);

	std::vector<double> Q(N);

	Dune::Timer watch;


	// Do sum initial samples

	for (int i = 0; i < N; i++){

		z.generate_random_field();

		Q[i] = model.getSample(level,z);
	}



	double time = watch.elapsed() / N;

	double V = computeVar(Q);

	double samplingError = std::sqrt(V/N);

	bool flag = false;

	if (samplingError > epsilon){flag = true;}

	int Nold;

	while(flag && N < maxSamples){

		Nold = N;
		N *= greedyFactor; N += 1;

		Q.resize(N);

		for (int i = Nold; i < N; i++){
			z.generate_random_field();
			Q[i] = model.getSample(level,z);
		}

		V = computeVar(Q);
		samplingError = std::sqrt(V/N);

		if (samplingError < epsilon){flag = false;}

	}

	double EQ = computeMean(Q);

	std::cout << "== Results == " << EQ << std::endl;
	std::cout << "E[Q] = " << EQ << std::endl;
	std::cout << "V[Q] = " << V << std::endl;
	std::cout << "Sampling Error = " << samplingError << std::endl;
	std::cout << "Average Time per Samples = " << time << std::endl;


}

	

