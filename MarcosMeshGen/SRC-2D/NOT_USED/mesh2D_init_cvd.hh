#ifndef _MESH_2D_INIT_HH
#define _MESH_2D_INIT_HH

#include<iostream>
using namespace std ;

#include<fstream>
#include<cassert>
#include<cstdlib>

#include<vector>
#include<cmath>
#include<string>

// general aux macros
#include "/Users/marco/Work/C++/Basics/my_macro.hh"

// mesh manager MFD applications
// (should be already called in main)
#include "./mesh2D.hh"

// mesh builder
#include "./mesh2D_builder.hh"

// functions for mesh building
#include "./mesh2D_greader.hh"
#include "./mesh2D_builtin.hh"
#include "./mesh2D_dualize.hh"
#include "./mesh2D_refiner.hh"

// curved faces 
#include "./mesh2D_cvd_builtin.hh"

// driver for functionalities
#include "./mesh2D_rpars.hh"

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
// RPARS is to be included before in main file !
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
// INITIALIZATION METHODS
// 
// flag = abs( mesh_flag )
//
// mesh_flag: 1..99:
// input from external data
// flag = 1    ==> read from input file mesh_name file1
// ----
// mesh_flag: Nxx, N=1..6
// internal construction
// flag = 1xx ==> quad mesh
// flag = 2xx ==> tria mesh
// flag = 3xx ==> Lshp mesh
// flag = 4xx ==> hexg mesh
// flag = 5xx ==> rotated triangles 
// flag = 6xx ==> double-quad grids (test case BLS-2)
// flag = 8xx ==> Octo mesh
// flag = 9xx ==> Octo mesh, non randomized
//
// mesh_flag < 0  ==> initialize dual mesh
//
//-------------------------------------------------------------------------------------------
// CONSTRUCTION OF PRIMAL MESH
// ilev>0 ==> internal refinement
void internal_mesh_generation( int flag, int ilev, mesh_2Dc & mesh, Read_Mesh2D_Pars & rpar ) {
  assert( 100<=flag && flag<=999 ) ;

  int    mf = flag-100*(flag/100)     ;
  int    nl = int( pow(2.,ilev) )     ; 
  int    nx = rpar.get_mesh_nx() * nl ; 
  int    ny = rpar.get_mesh_ny() * nl ; 
  double wx = rpar.get_mesh_wx()      ; 
  double wy = rpar.get_mesh_wy()      ; 

  PRT(flag/100) ;
  PRT(mf) ;

  switch( flag/100 ) {
  case 1: build_new_quad_mesh     ( mesh, nx, ny, mf, wx, wy ) ; break ; // 1xx --> quads
  case 2: build_tria_mesh         ( mesh, nx, ny, mf, wx, wy ) ; break ; // 2xx --> triangles
  case 3: build_lshape_mesh       ( mesh, nx, ny, mf, wx, wy ) ; break ; // 3xx --> L-shape
    //case 4: build_hexa_mesh         ( mesh, nx                 ) ; break ; // 400 --> hexagons
  case 5: build_rotated_tria_mesh ( mesh, nx                 ) ; break ; // 500 --> rot-trias
  case 6: build_double_quad_mesh  ( mesh, nx, ny, mf, wx, wy ) ; break ; // 6xx --> double-quads
  case 7: build_curved_face_mesh  ( mesh, nx, ny, mf, wx, wy, ilev ) ; break ; // 7xx --> double-quads
  case 8: build_NonConvexOcto_mesh( mesh, nx, ny, mf, wx, wy ) ; break ; // 8xx --> octagons
  case 9: build_NonConvexNonRandomOcto_mesh( mesh, nx, ny, mf, wx, wy ) ; break ; // 9xx --> octagons, non randomized
  default: assert(false) ;
  } 
}
void read_input_mesh( mesh_2Dc & mesh, Read_Mesh2D_Pars & rpar ) {
  int default_offset = FORTRAN_OFFSET ;// OFFSET ;
  string file_name = rpar.get_mesh_name() ;
  mesh2D_reader_GeneralFormat input_mesh( mesh, file_name, default_offset ) ;
  input_mesh.read_and_build() ;
}
// this mesh dualization works only if primal cells are convex!
void dual_mesh_generation( mesh_2Dc & mesh, 
			   mesh_2Dc & dual_mesh, Read_Mesh2D_Pars & rpar ) {
  bool tria_based_mesh = true ;
  bool quad_based_mesh = true ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    tria_based_mesh &= mesh.n_regn_face(iR)==3 ;
    quad_based_mesh &= mesh.n_regn_face(iR)==4 ;
  }
  if ( tria_based_mesh || quad_based_mesh ) {
    build_dual_mesh( mesh, dual_mesh ) ;
  } else {
    cout << "ERROR MESSAGE: Dualization is possible only for meshes of triangles or quads" << endl ;
    assert(false) ;
  }
}
// regular refinement for triangles and quads
void build_regular_refinement( mesh_2Dc & mesh, mesh_2Dc & split_mesh,
			       Read_Mesh2D_Pars & rpar ) {
  bool tria_based_mesh = true ;
  bool quad_based_mesh = true ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    tria_based_mesh &= mesh.n_regn_face(iR)==3 ;
    quad_based_mesh &= mesh.n_regn_face(iR)==4 ;
  }
  if      ( tria_based_mesh ) { build_split_triangular_mesh( mesh, split_mesh ) ; }
  else if ( quad_based_mesh ) { assert(false) ; } //{ build_split_polygonal_mesh ( mesh, split_mesh ) ; }
  else                        { assert(false) ; } 
}
/*
  -- sign of mesh_flag
  mesh_flag>0 ==> works on primal mesh
  mesh_flag<0 ==> works on dual   mesh

  -- inpt_flag = abs(mesh_flag)
  inpt_flag<100                    ==> first mesh is input from disk files
  100<=inpt_flag && inpt_flag<=600 ==> all meshes are internally generated

  -- refinement strategy:
  -  ====================
  -  if first mesh is from external input
  -     ==> only regular refinement (independent of refn_flag)
  
  -  if first mesh is from internal routines 
  -     ==> choose between:
  -         (i)  homothetic refinement
  -         (ii) new internal generation

  -- homothetic refinement automatically selects between
  -         (i)   nested triangles if triangle-based mesh    (theta=1)
  -         (ii)  nested quads     if quad-based     mesh    (theta=0)
  -         (iii) rotated homothetic refn if generic polygon (theta in (0,1))
  -
  - theta controls the size of the new cell inside each region
  - better values of theta are in the range [0.45,0.55]

  -- dualization (for mesh_flag<0):
  -  first mesh is always dualized
  -  refined meshes are dualized only if primal mesh is internally generated

*/
void init_mesh( int ilev, mesh_2Dc & mesh, Read_Mesh2D_Pars & rpar ) {
  assert( ilev>=0 ) ;
  // flags
  int mesh_flag = rpar.get_mesh_flag() ;
  int refn_flag = rpar.get_refn_flag() ;
  int inpt_flag = abs(mesh_flag) ;
  assert( 0<inpt_flag && inpt_flag<=999 ) ;
  
  if ( inpt_flag<100 ) {
    // mesh generation from EXTERNAL data on disk files
    if ( ilev==0 ) { read_input_mesh( mesh, rpar ) ; }
    else           { 
      read_input_mesh( mesh, rpar ) ;
      for ( int i=0; i<ilev; ++i ) {
	LINE(---) ;
	MSGF(endl<<"Compute refinement level --> "<<i<<endl) ;
	build_regular_refinement ( mesh, mesh, rpar ) ; 
	MSGF(endl<<"end of refinement level --> "<<i<<endl) ;
      }
    }
  } else if ( 100<=inpt_flag && inpt_flag<=999 ) { 
    // mesh generation from INTERNAL routines
    if      ( ilev>0 && refn_flag==0 || ilev==0 ) { internal_mesh_generation( inpt_flag, ilev, mesh, rpar ) ; }
    else if ( ilev>0 && refn_flag==1            ) { build_regular_refinement( mesh, mesh, rpar ) ;            } 
  } else {
    assert(false) ;
  }
  // dualization if required (and possible) for triangle- and quad-based meshes
  if ( mesh_flag<0 ) {
    dual_mesh_generation( mesh, mesh, rpar ) ;
  }
} // end of void init_mesh( ...

// only for generating dual mesh from homothetically refined mesh
void init_mesh( int ilev, mesh_2Dc & mesh, mesh_2Dc & dual_mesh, Read_Mesh2D_Pars & rpar ) {
  assert( ilev>=0 ) ;
  // flags ---
  int mesh_flag = rpar.get_mesh_flag() ;
  int refn_flag = rpar.get_refn_flag() ;
  int inpt_flag = abs(mesh_flag) ;
  assert( 0<inpt_flag && inpt_flag<=999 ) ;

  if ( inpt_flag<100 ) {
    // mesh generation from EXTERNAL data on disk files
    if ( ilev==0 ) { read_input_mesh( mesh, rpar ) ; }
    else           { 
      read_input_mesh( mesh, rpar ) ;
      for ( int i=0; i<ilev; ++i ) {
	LINE(---) ;
	MSGF(endl<<"Compute refinement level --> "<<i<<endl) ;
	build_regular_refinement ( mesh, mesh, rpar ) ; 
	MSGF(endl<<"end of refinement level --> "<<i<<endl) ;
      }
    }
  } else if ( 100<=inpt_flag && inpt_flag<=999 ) { 
    // mesh generation from INTERNAL routines
    if ( ilev==0 || ilev>0  && refn_flag==0 ) {
      internal_mesh_generation( inpt_flag, ilev, mesh, rpar ) ;
    } else if ( ilev>0 && refn_flag==1 ) {
      build_regular_refinement( mesh, mesh, rpar ) ;
    }
  } else {
    assert(false) ;
  }
  // dualization if required (and possible) for triangle-based meshes
  if ( mesh_flag<0 ) {
    dual_mesh_generation( mesh, dual_mesh, rpar ) ;
  }
} // end of void init_mesh( ...

#endif // end of _MESH_2D_INIT_HH
