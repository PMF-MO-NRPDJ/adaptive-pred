#pragma once

#include <cmath>
#include <iostream>
#include <string>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/gmshreader.hh>
//#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/pdelab/adaptivity/adaptivity.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

//#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
//#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include "operator.hh"
#include "estimator.hh"
#include "bctype.hh"


template <typename Grid>
void driver(Grid &grid)
{
  using GV = typename Grid::LeafGridView;
  GV gv = grid.leafGridView();
  const int dim = GV::dimension;
  // P1 polinomi na simpleksima
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV, double, double, 1>;
  FEM fem(gv);
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  DirichletBdry bctype; // identifikacija Dirichletove granice
  using CC = typename GFS::template ConstraintsContainer<double>::Type;
  CC cc;
  Dune::PDELab::constraints(bctype, gfs, cc);

  // Vektor stupnjeva slobode rješenja
  using U = Dune::PDELab::Backend::Vector<GFS, double>;
  U u(gfs);
  BCExtension<GV> bcext(gv);
  // Izračunaj z iz rubnog uvjeta
  Dune::PDELab::interpolate(bcext, gfs, u);

  // Petlja unutar koje adaptiramo mrežu.
  double tol = 0.05;   // Tolerancija za procjenitelj
  double alpha = 0.8;  // 0 <= alpha <= 1
  double beta = 0.0;   // 0<= beta <= 1
  int no_steps = 5;    // max broj koraka profinjenja

  for (int i = 0; i < no_steps; ++i)
  {
    std::string iter = std::to_string(i);
    std::cout << "===== Iteracija no. : " << iter
              << "  Najviše profinjenje mreže: " << grid.maxLevel() << std::endl;
    std::cout << "      Br vezanih stupnjeva slobode = "  << cc.size() << " od " << gfs.globalSize()
              << std::endl;

    // Najsigurnije je imati definicije prostora unutar petlje modifikacije mreže.
    // Lokalni operator
    using LOP = LocalOperator<DirichletBdry,FEM>;
    LOP lop(bctype);

    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe( 9 ); // Za P1 elemente 3^dim
    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
    GO go(gfs, cc, gfs, cc, lop, mbe);

    using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO>;
    LS ls(100, 0);

    using  SLP =Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U>;
    // Nađi aproksimativno rješenje.
    const double redukcija = 1e-10;
    const double min_defect = 1e-99;
    int verbosity = 0;
    SLP slp(go, ls, u, redukcija, min_defect, verbosity);
    slp.apply();

    // Konstrukcija procjenitelja greške.
    // dopuni ....
    //
    // Formiranje indikatora greške.
    // dopuni ....   
    auto estimated_error = sqrt(indicator.one_norm()); // L^1 norma jer su u indicator kvadrati
    std::cout << "      L2 norma greške je procijenjena na " << estimated_error << std::endl;

    // VTK ispis
    // Rješenje
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    using UDGF = Dune::PDELab::DiscreteGridFunction<GFS, U>;
    UDGF udgf(gfs, u);
    using VTKF = Dune::PDELab::VTKGridFunctionAdapter<UDGF>;
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(udgf, "fem_sol")));

    // Indikator greške
    // dopuni ...

    std::string output = "adaptive";
    vtkwriter.write(output + iter, Dune::VTK::ascii);

    // Da li je greška dovoljno mala ?
    if (estimated_error <= tol)
      break;  // prekini profinjavanje
    if (i == no_steps - 1)
      break; // u zadnjem koraku preskoči adaptaciju mreže

    // Označiti elemente za profinjenje i adaptiraj mrežu
    // dopuniti     ....

    // ponovo izračunaj Dirichletova ograničenja
    Dune::PDELab::constraints(bctype, gfs, cc);
    // korektne rubne uvjete upiši u novi vektor
    U unew(gfs);
    Dune::PDELab::interpolate(bcext, gfs, unew);
    // kopiraj Dirichletove rubne uvjete u interpolirani vektor rješenja
    Dune::PDELab::copy_constrained_dofs(cc, unew, u);
  }
}

