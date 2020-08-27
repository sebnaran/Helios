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

void postscript_output( mesh_2Dv & mesh ) {

  // set file name
  string psout_file = string("./MESH-PS/mesh_2D") ;
  
  // print postscript picture of mesh
  PS_Driver psd(psout_file) ;
  mesh2D_psout msh_psout( mesh, psd ) ;

  // just to change offset
  msh_psout.set_offset( OFFSET ) ;
  
  // draw the mesh
  bool draw_numbers = false ;
  string mesh_name = mesh.get_mesh_name() ;
  if ( mesh_name==string("") ) { mesh_name=string("Mesh") ; }
  msh_psout . draw_mesh(mesh_name,draw_numbers) ;

  // bullets?
  bool draw_bullets = false ;  
  if ( draw_bullets ) {
    msh_psout . draw_vrtx_bullet( 3 ) ; 
  }
}

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

  // print all log files
  if ( true ) {
    mesh2Dv_printer mesh_printer(mesh) ;
    mesh_printer.print_all_datasets() ;
    mesh_printer.print_all_regions () ;
    mesh_printer.print_all_faces   () ;
    mesh_printer.print_all_vertices() ;
    mesh_printer.print_additional_data() ;
  }
}
