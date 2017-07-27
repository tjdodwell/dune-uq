// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

//#include "defaultimp.hh"
//#include "pattern.hh"
#include "flags.hh"
//#include "idefault.hh"

namespace Dune {
    namespace PDELab {

        template <int dim>
        class AverageInSubDomain : public NumericalJacobianApplyVolume<AverageInSubDomain<dim> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<AverageInSubDomain<dim> >
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };

            AverageInSubDomain (Dune::FieldVector<double,dim>& y_, double radius_,int intorder_=5)
            : y(y_),radius(radius_),intorder(intorder_)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {

                // domain and range field type
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename R::value_type RF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::JacobianType JacobianType;
                typedef typename LFSU::Traits::SizeType size_type;

                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeType RangeType;

                // select quadrature rule
                GeometryType gt = eg.geometry().type();
                const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

                // Loop over quadrature points
                for (const auto& ip : rule)
                {
                  // evaluate basis functions
                  std::vector<RangeType> phi(lfsu.size());
                  lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);


                  // evaluate u at integation points
                   double u = 0.0;
                   Dune::FieldVector<double,dim> xip = eg.geometry().global(ip.position());

                   for (int i = 0; i < phi.size(); i++){
                       u += phi[i] * x(lfsu,i);
                   }

                  double d = 0.0;
                  for (int j = 0; j < dim; j++){ d += (xip[j] - y[j]) * (xip[j] - y[j]); }
                  d = std::sqrt(d);


                  double factor = 0.0;
                  if (d < radius){ factor = ip.weight() * eg.geometry().integrationElement(ip.position());}

                    for (int i = 0; i < lfsu.size(); i++){
                          r.accumulate(lfsv,i,u * phi[i] * factor);
                    }

                } // end for each quadrature point

            } // end alpha_volume


        private:
            int intorder;
            Dune::FieldVector<double,dim>& y;
            double radius;

        };


}
}
