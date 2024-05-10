#pragma once
#include<vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "coefficients.hh"

/** Lokalni operator za eliptičku zadaću
 *
 *   -\Delta u(x) + a u(x) = f(x) x\in\Omega,
 *                     u(x) = g(x) x\in\partial\Omega
 *
 */
template<typename DirichletBdry, typename FEM>
class LocalOperator :
  public Dune::PDELab::NumericalJacobianVolume<LocalOperator<DirichletBdry,FEM> >,
  public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<DirichletBdry,FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  DirichletBdry const & dirBdry;

public:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  LocalOperator (DirichletBdry const & dirBdry_) : dirBdry(dirBdry_) {}

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim = EG::Entity::dimension;
    using Gradient = Dune::FieldVector<double,dim>;

    // kvadraturna formula
    auto geo = eg.geometry();
    const int order = 2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    for (const auto& ip : rule)
      {
        auto xi = ip.position();
        auto glob_xi = eg.geometry().global(xi);
        auto& phihat = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis());

        // evaluate u
        double  u=0.0;
        for (size_t i=0; i<lfsu.size(); i++)
          u += x(lfsu,i)*phihat[i];

        auto& gradphihat = cache.evaluateJacobian(xi, lfsu.finiteElement().localBasis());
        const auto jac = geo.jacobianInverseTransposed(xi);
        std::vector<Gradient> gradphi(lfsu.size());

        for (size_t i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        // grad u
        Dune::FieldVector<double,dim> gradu(0.0);
        for (size_t i=0; i<lfsu.size(); i++)
          gradu.axpy(x(lfsu,i),gradphi[i]);

        auto factor = ip.weight()*geo.integrationElement(xi);
        auto a = reactCoeff(glob_xi);
        auto f= RHS(eg.geometry().global(xi));

        for (size_t i=0; i<lfsu.size(); i++)
          r.accumulate(lfsu,i,( gradu*gradphi[i] + a*u*phihat[i] - f*phihat[i]) * factor);
      }
  }
};
