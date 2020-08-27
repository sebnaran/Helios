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
  mesh_2Dv & mesh_0 = rpar.get_mesh_flag()>0 ? primal_mesh : dual_mesh ;

  // build the base mesh on the square domain
  int nlev  = rpar.get_mesh_nlev() ;
  init_mesh( nlev, primal_mesh, dual_mesh, rpar ) ;

  // set the domain size
  double Lx = 1 ;
  double Ly = 1 ;

  // rescale the domain in a [0,Lx]x[0,Ly] rectangle
  double xmin = 0 ;
  double ymin = 0 ;
  double xmax = 1 ;
  double ymax = 1 ;
  mesh_0.change_bbox( xmin, ymin, xmax, ymax ) ;
  post_proc( nlev, mesh_0, rpar ) ;

  mesh2Dv_writer mesh_writer(mesh_0) ;
  mesh_writer.write_mesh_Durham_format("locrefs_quads_0") ;

  // refined meshes
  mesh_2Dv mesh_1, mesh_2 ;
  
  { // first refinement: mesh_0 --> mesh_1
    vector<int> flagged_cells_list ;
    double xb = 0.1*Lx ;
    double xe = 0.9*Lx ;
    double yb = 0.1*Ly ;
    double ye = 0.9*Ly ;
    int nR = mesh_0.n_region() ;
    for ( int iR=0; iR<nR; ++iR ) {
      double xR = mesh_0.coords_R(iR,0) ;
      double yR = mesh_0.coords_R(iR,1) ;
      bool do_refinement = ( xb<=xR && xR<=xe ) && ( yb<=yR && yR<=ye ) ;
      if ( do_refinement ) {
	flagged_cells_list.push_back(iR) ;
      }
    }

    int nrx = 1 ; // refinement rate in x
    int nry = 1 ; // refinement rate in y

    LocalMeshRefinement loc_msh_ref(mesh_0) ;
    loc_msh_ref.setup( flagged_cells_list, nrx, nry ) ;

    loc_msh_ref.build_local_mesh_refinement( mesh_1 ) ;
    post_proc( nlev+1, mesh_1, rpar ) ;

    mesh2Dv_writer mesh_writer(mesh_1) ;
    mesh_writer.write_mesh_Durham_format("locrefs_quads_1") ;
  }

  { // second refinement: mesh_1 --> mesh_2
    vector<int> flagged_cells_list ;
    double xb = 0.1*Lx ;
    double xe = 0.9*Lx ;
    double yb = 0.3*Ly ;
    double ye = 0.7*Ly ;
    int nR = mesh_1.n_region() ;
    for ( int iR=0; iR<nR; ++iR ) {
      double xR = mesh_1.coords_R(iR,0) ;
      double yR = mesh_1.coords_R(iR,1) ;
      bool do_refinement = ( xb<=xR && xR<=xe ) && ( yb<=yR && yR<=ye ) ;
      if ( do_refinement ) {
	flagged_cells_list.push_back(iR) ;
      }
    }

    int nrx = 1 ; // refinement rate in x
    int nry = 1 ; // refinement rate in y

    LocalMeshRefinement loc_msh_ref(mesh_1) ;
    loc_msh_ref.setup( flagged_cells_list, nrx, nry ) ;

    loc_msh_ref.build_local_mesh_refinement( mesh_2 ) ;
    post_proc( nlev+2, mesh_2, rpar ) ;

    mesh2Dv_writer mesh_writer(mesh_2) ;
    mesh_writer.write_mesh_Durham_format("locrefs_quads_2") ;
  }
  
  // print all log files
  if ( false ) {
    mesh2Dv_printer mesh_printer(mesh_2) ;
    mesh_printer.print_all_datasets() ;
    mesh_printer.print_all_regions () ;
    mesh_printer.print_all_faces   () ;
    mesh_printer.print_all_vertices() ;
    mesh_printer.print_additional_data() ;
  }

}
