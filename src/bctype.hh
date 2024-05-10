#pragma once

#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/common/function.hh>


#include "coefficients.hh"

class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
      return true;
  }

};

template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >,
                                      BCExtension<GV> > 
{
  GV const & gv;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >;

  BCExtension (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
      auto x = e.geometry().global(xlocal);     
      y =  exact(x);
  }

  const GV& getGridView () const {return gv;}
};




