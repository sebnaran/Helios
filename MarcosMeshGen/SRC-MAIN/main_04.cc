#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <complex>
using namespace std ;

// BASIC FILES
#include "../SRC-MISC/my_macro.hh"

// MFD RUN PARS
#include "../SRC-MISC/read_parameters.hh"

// 2D MESH
const int FORTRAN_OFFSET = 1 ;
const int CPP_OFFSET     = 0 ;
const int OFFSET = CPP_OFFSET ;
#include "../SRC-2D/mesh2De.hh"
#include "../SRC-2D/mesh2D_init.hh"
#include "../SRC-2D/mesh2D_post.hh"
#include "../SRC-2D/mesh2D_writer.hh"
#include "../SRC-2D/mesh2D_localrefiner.hh"

int main() {

  // read input parameters
  RPars rpar("data.inp") ;
  rpar.read_file() ;
  rpar.print_parameters() ;
  
  // instantiate [primal|dual] (empty) mesh variable
  mesh_2Dv primal_mesh, dual_mesh ;
  mesh_2Dv & mesh = rpar.get_mesh_flag()>0 ? primal_mesh : dual_mesh ;

  int nlev  = rpar.get_mesh_nlev() ;
  init_mesh( nlev, primal_mesh, dual_mesh, rpar ) ;
  post_proc( nlev, mesh, rpar ) ;

  vector<int> flagged_cells_list ;
#if 0
  flagged_cells_list.push_back(0) ;
  flagged_cells_list.push_back(2) ;
  flagged_cells_list.push_back(1) ;
#endif

  int nrx = 3 ;
  int nry = 2 ;

  LocalMeshRefinement loc_msh_ref(mesh) ;
  loc_msh_ref.setup( flagged_cells_list, nrx, nry ) ;

  mesh_2Dv new_mesh ;
  loc_msh_ref.build_local_mesh_refinement( new_mesh ) ;
  post_proc( nlev+1, new_mesh, rpar ) ;

  // print all log files
  if ( true ) {
    mesh2Dv_printer mesh_printer(new_mesh) ;
    mesh_printer.print_all_datasets() ;
    mesh_printer.print_all_regions () ;
    mesh_printer.print_all_faces   () ;
    mesh_printer.print_all_vertices() ;
    mesh_printer.print_additional_data() ;
  }

}
