#include <random>

#include "../UQ/KLFunctions.h"


class RandomField
{

public:

	RandomField(){

		// KL parameters used for 2D tests
    	sigKL_ = config.get<double>("RandomField.sig",2.5);
    	ellKL_ = config.get<double>("RandomField.ell",1./3.);
    	numKLmodes_ = config.get<int>("RandomField.KLmodes",100);

        const int dim = 2;

        xi.resize(numKLmodes_);

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

    double inline evaluateScalar(const Dune::FieldVector<double,2> x) const {	

        double k = 0.0;
        for (int j = 0; j < numKLmodes_; j++){
            k += std::sqrt(lam[j]) * evalPhi(x[0],mode_id_i[j],1.0) * evalPhi(x[1],mode_id_j[j],1.0) * xi[j];
        }
        k = std::exp(k);

        return k;	

    }

	void inline generate_random_field(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> randn(0.0,sigKL_);
        std::fill(xi.begin(), xi.end(), 0.0);
        for (int j = 0; j < numKLmodes_; j++){
            xi[j] = randn(gen);
        }       
    } // end user_random_field

    double inline getValue(int i){
    	return xi[i];
    }


private:
	
	int numKLmodes_;
    bool isTest = false;
    double sigKL_, ellKL_;
    std::vector<double> xi;
    std::vector<double> freq, lam1D, lam, param_;
    std::vector<int> mode_id_i, mode_id_j;
    std::vector<Dune::FieldMatrix<double,2,2>> Kij;


}; // end of RandomField class