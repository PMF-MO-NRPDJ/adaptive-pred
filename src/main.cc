#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/uggrid.hh>

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include "driver.hh"

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 2;
  std::string filename = "ldomain.msh";

  using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
  std::unique_ptr<Grid> gridp = Dune::GmshReader<Grid>::read(filename);

  driver(*gridp);

  return 0;
}
