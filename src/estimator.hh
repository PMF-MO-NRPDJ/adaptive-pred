#pragma once
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

//#include<dune/pdelab/common/referenceelements.hh> 2.6.0
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "coefficients.hh"


/** Rezidualni procjenitelj za eliptičku zadaću:
 *
 *   -\Delta u(x) + a u(x) = f(x) x\in\Omega,
 *                     u(x)  = g(x) x\in\partial\Omega
 *
 *
 * Za izračunati procjenitelj treba pozvati residual() na mrežnom
 * operatoru. Lokalno se računa \eta_K^2 na svakom elementu K.
 *
 * Pretpostavke i ograničenja:
 * - Uzima se da je  LFSU jednak P_k/Q_k te da je
 *   LFSV jednak P_0.
 * - Derivacije drugog reda se zanemaruju.
 *
 */
template<typename BCType, typename FEM>
class Estimator : public Dune::PDELab::LocalOperatorDefaultFlags
{
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  BCType const & bctype;

  // dijametar ćelije
  template<class GEO>
  double diameter(const GEO& geo) const
  {
    double hmax = -1.0E00;
    for (int i=0; i<geo.corners(); i++)
      {
        auto xi = geo.corner(i);
        for (int j=i+1; j<geo.corners(); j++)
          {
            auto xj = geo.corner(j);
            xj -= xi;
            hmax = std::max(hmax,xj.two_norm());
          }
      }
    return hmax;
  }

public:
  enum { doPatternVolume = false };
  enum { doPatternSkeleton = false };
  enum { doAlphaVolume  = true };
  enum { doAlphaSkeleton  = true };

  Estimator(BCType const & bctype_) : bctype(bctype_) {}

  // volumni dio indikatora
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {

    const int dim = EG::Geometry::coorddimension;
    auto geo = eg.geometry();
    const int order = 2;
    auto & rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),order);
    // jednostavnija mogućnost
    //  auto rule = Dune::PDELab::quadratureRule(geo,order);

    double sum=0.0;
    for (const auto& ip : rule)
      {
        auto xi = ip.position();
        auto glob_xi = eg.geometry().global(xi);
        auto& phihat = cache.evaluateFunction(xi,lfsu.finiteElement().localBasis());

        // nastaviti    
      }

  // Dio indikatora po unutarnjim stranicama elemenata.
  // Svaka se stranica obilazi samo jednom.
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
      const LFSU & lfsu_i, const X & x_i, const LFSV & lfsv_i,
      const LFSU & lfsu_o, const X & x_o, const LFSV & lfsv_o,
      R& r_i, R& r_o) const
  {
    const int dim = IG::Geometry::coorddimension;

    auto globalgeo = ig.geometry();
    const int order = 2;
    auto & rule = Dune::QuadratureRules<double,dim-1>::rule(globalgeo.type(),order);
    // jednostavnija verzija
    //auto rule = Dune::PDELab::quadratureRule(globalgeo,order);
    double sum=0.0;
    for(const auto& ip : rule)
    {
         // nastaviti        
    }

    // akumulacija indikatora
    auto h_T = diameter(globalgeo);
    r_i.accumulate(lfsv_i, 0, 0.5*h_T * sum);
    r_o.accumulate(lfsv_o, 0, 0.5*h_T * sum);
  }

 };
