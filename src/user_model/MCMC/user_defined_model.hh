#ifndef user_model_h
#define user_model_h


#include "diffusion.hh"
#include "AverageInSubDomain.hh"

// ------------------------------------------------------------------------
//             Dirichlet Boundary Conditions Class
// ------------------------------------------------------------------------
//! \brief Parameter class selecting boundary conditions
class BCTypeParam
: public Dune::PDELab::DirichletConstraintsParameters
{
public:
    //! Test whether boundary is Dirichlet-constrained
    template<typename I>
    bool isDirichlet(const I & intersection,
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
    {
        Dune::FieldVector<typename I::ctype, I::dimension>
        xg = intersection.geometry().global( coord );

        return true;  // Dirichlet b.c. on all other boundaries
    }

};

template<int dim, typename GRID>
class MODEL{

	public:

		MODEL(GRID& grid) : grid_(grid){

            sigF = config.get<double>("data.sigF",1e-2);


        }

		double inline getSample(int l, RandomField& z, bool proposal = true){

			using Dune::PDELab::Backend::native;

			Dune::Timer watch;

			typedef typename GRID::LevelGridView LGV;
        	LGV gv = grid_.levelGridView(l);

        	typedef double RF;
        	typedef typename LGV::Grid::ctype Coord;

        	const int dofel = 4;

	        typedef Dune::PDELab::QkLocalFiniteElementMap<LGV,Coord,RF,1> FEM;
	        FEM fem(gv);

	        typedef Dune::PDELab::ConformingDirichletConstraints CON;

	        typedef Dune::PDELab::istl::VectorBackend<> VectorBackend;

	        typedef Dune::PDELab::GridFunctionSpace<LGV, FEM, CON, VectorBackend> GFS;
	        GFS gfs(gv,fem); gfs.name("pressure");

	        // <<<3>>> assemble constraints on this space
            BCTypeParam bctype;

            typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
                CC cc;

            Dune::PDELab::constraints( bctype, gfs, cc );

            typedef typename Dune::PDELab::Backend::impl::BackendVectorSelector<GFS,RF>::Type U;
                U p(gfs,0.0);

	        typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
	        MBE mbe(9); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

            int intOrder = config.get<int>("solver.intOrder",2);

	        typedef Dune::PDELab::diffuse<RandomField,LGV> LOP;
	        LOP lop(gv,z,config.get<double>("PDE.f",1.0),intOrder,proposal);

	        typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC> GO;
	        GO go(gfs,cc,gfs,cc,lop,mbe);

	        // Select a linear solver backend
	        typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
	        LS ls(5000,0);

	        // Select linear problem solver
	        typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
	        SLP slp(go,ls,p,1e-10,1e-99,0);

	        slp.apply(); // Compute solution to test problem


	        if(config.get<bool>("vtk.solution",false)){
	          Dune::VTKWriter<LGV> vtkwriter(gv,Dune::VTK::conforming);
	          Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,p);
	          vtkwriter.write("solution",Dune::VTK::appendedraw);
	        }
	    
	        // Compute Quantity of Interest

            Dune::FieldVector<double,2> y(0.0); 

            y[0] = config.get<double>("QoI.center_x",0.75);
            y[1] = config.get<double>("QoI.center_y",0.75);

            double radius = config.get<double>("QoI.radius",0.15);;

            //  Construct Linear Operator on FEM Space
            typedef Dune::PDELab::AverageInSubDomain<2> QLOP;
                QLOP Qlop(y,radius);

            typedef Dune::PDELab::GridOperator<GFS,GFS,QLOP,MBE,double,double,double,CC,CC> QGO;
                QGO qgo(gfs,cc,gfs,cc,Qlop,mbe);
            U qq(gfs,0.0);

            qgo.residual(p,qq);

            double Q = 0.0;
            for (int ii = 0; ii < native(qq).size(); ii++){
                Q += native(qq)[ii];
            }

            Q /=  M_PI * radius * radius;

            // Compute Likelihood

            // Compute Data first

            Dune::FieldVector<double,3> Q_data(0.0);

            std::vector<int> nodal_data = config.get<std::vector<int> >("data.nodes",{33,58,82});

            double likelihood = 0.0;
            for (int i = 0; i < numdata; i++){
                likelihood += (native(p)[nodal_data[i]] - data[i]) * (native(p)[nodal_data[i]] - data[i]) 
            }

            likelihood = std::exp(-likelihood/sigF);

            z.set_likelihood(likelihood,proposal);
	       
            return Q;

		}


	private:

        GRID& grid_;
        double sigF;


};





#endif /* user_model_h */



