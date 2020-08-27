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
const int OFFSET = FORTRAN_OFFSET ;
#include "../SRC-2D/mesh2De.hh"
#include "../SRC-2D/mesh2D_init.hh"
#include "../SRC-2D/mesh2D_post.hh"
#include "../SRC-2D/mesh2D_writer.hh"

// only from input file using Durham format
// can use different construction

int main() {

  // read input parameters
  RPars rpar("data.inp") ;
  rpar.read_file() ;
  rpar.print_parameters() ;

  string file_name = rpar.get_mesh_name() ;
  int default_offset = OFFSET ; // FORTRAN_OFFSET ;// OFFSET ;
  //mesh2D_reader_GeneralFormat input_mesh( mesh, file_name, default_offset ) ;

  // instantiate mesh
  mesh_2Dv mesh ;

  int mesh_flag = rpar.get_mesh_flag() ;
  mesh2D_reader_DurhamFormat input_mesh( mesh, file_name, default_offset ) ;
  input_mesh.read_and_build( mesh_flag ) ;

  // post-script output
  int nlev = rpar.get_mesh_nlev() ;
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
