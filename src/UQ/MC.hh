#include "../stats/common.hh"

/* MC.hh - Simple Monte Carlos Algorithms */
template <class MODEL, class RandomField, class Samples>
void inline MC(MODEL& model, RandomField& z,Samples& samples, int level = 0){

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


	// Do some initial samples
	for (int i = 0; i < N; i++){

		watch.reset();

		z.generate_random_field();

		Q[i] = model.getSample(level,z);

		samples.add_Sample(Q[i],z, watch.elapsed());
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

	

