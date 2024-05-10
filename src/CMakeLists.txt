add_executable(main main.cc driver.hh
               estimator.hh
               operator.hh
               bctype.hh
               coefficients.hh
           )


dune_symlink_to_source_files(FILES ldomain.msh)
