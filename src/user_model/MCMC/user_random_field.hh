#include <random>

#include "../../UQ/KLFunctions.h"


class RandomField
{

public:

	RandomField(){

		// KL parameters used for 2D tests
    	sigKL_ = config.get<double>("RandomField.sig",2.5);
    	ellKL_ = config.get<double>("RandomField.ell",1./3.);
    	numKLmodes_ = config.get<int>("RandomField.KLmodes",100);

        sigPCN_ = config.get<double>("Proposal.sigPCN",1.0);
        beta = config.get<double>("Proposal.beta",0.1);

        const int dim = 2;

        xi.resize(numKLmodes_);
        xi_p.resize(numKLmodes_);

        // Define and Resize vector for KL modes definition
        int N = std::pow(numKLmodes_,1.0/dim) + 1;
        freq.resize(N); lam1D.resize(N); lam.resize(numKLmodes_);
        mode_id_i.resize(numKLmodes_);
        mode_id_j.resize(numKLmodes_);

        rootFinder(N,ellKL_,freq);
        evaluate_eigenValues(ellKL_,sigKL_,lam1D,freq);
        construct_2D_eigenValues(lam1D, lam, mode_id_i, mode_id_j);

	};

	double inline evalPhi(double x,int i, double L) const {
        double phi = 0.0;
        double omega = freq[i];
        x -= L;
        double tmp = sqrt(L + std::sin(2 * omega * L) / (2 * omega));
        if (i % 2 == 0){ phi = std::cos(omega * x) / tmp;}
        else { phi = std::sin(omega * x) / tmp; }
        return phi;
    }

    double inline evaluateScalar(const Dune::FieldVector<double,2> x, bool prop = false) const {	

        double k = 0.0;

        if (prop){
            for (int j = 0; j < numKLmodes_; j++){
                k += std::sqrt(lam[j]) * evalPhi(x[0],mode_id_i[j],1.0) * evalPhi(x[1],mode_id_j[j],1.0) * xi_p[j];
            }
        }
        else{
            for (int j = 0; j < numKLmodes_; j++){
                k += std::sqrt(lam[j]) * evalPhi(x[0],mode_id_i[j],1.0) * evalPhi(x[1],mode_id_j[j],1.0) * xi[j];
            }
        }
        k = std::exp(k);

        return k;	

    }

	void inline generate_random_field(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> randn(0.0,sigKL_);
        
        std::fill(xi.begin(), xi.end(), 0.0);

        prior_density = 0.0;

        for (int j = 0; j < numKLmodes_; j++){
            xi[j] = randn(gen);
            prior_density -= xi[j] * xi[j];
        }    

        double factor = std::pow((2. * M_PI), 0.5 * numKLmodes_);

        prior_density = std::exp(-0.5 * prior_density) / factor;



    } // end user_random_field

    void inline generate_proposal(){

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> randn(0.0,sigPCN_);

        prior_density_prop = 0.0;

        for (int j = 0; j < numKLmodes_; j++){
            xi_p[j] = std::sqrt(1 - beta * beta) * xi_p[j] + beta * randn(gen);
            prior_density -= xi_p[j] * xi_p[j];
        }

        double factor = std::pow((2. * M_PI), 0.5 * numKLmodes_);

        prior_density_prop = std::exp(-0.5 * prior_density_prop) / factor;


    }

    void inline set_likelihood(double l, bool prop = false){
        if (prop){ likelihood_p = l; }
        else {likelihood = l; }
    }

    bool inline accept_proposal(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> rand(0.0,1.0);
        double ratio = (prior_density_prop * likelihood_p) / (prior_density * likelihood);
        bool accept = false;
        if (rand(gen) < ratio){ // accept proposal
            for (int j = 0; j < numKLmodes_; j++){ xi[j] = xi_p[j]; }
            likelihood = likelihood_p;
            prior_density = prior_density_prop;
            accept = true;
        }
        return accept;
    }

    double inline getValue(int i){
    	return xi[i];
    }


private:
	
	int numKLmodes_;
    bool isTest = false;
    double likelihood, likelihood_p, prior_density, prior_density_prop;
    double sigKL_, ellKL_, sigPCN_, beta;
    std::vector<double> xi, xi_p;
    std::vector<double> freq, lam1D, lam, param_;
    std::vector<int> mode_id_i, mode_id_j;
    std::vector<Dune::FieldMatrix<double,2,2>> Kij;


}; // end of RandomField class