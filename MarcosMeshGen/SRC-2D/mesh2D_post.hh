#ifndef _MESH_2D_POST_HH
#define _MESH_2D_POST_HH

#include<iostream>
using namespace std ;

#include<fstream>
#include<cassert>
#include<cstdlib>

#include<vector>
#include<cmath>
#include<string>

// general aux macros
#include "../SRC-MISC/my_macro.hh"

// mesh manager for MFD applications
// (should be already called in main)
#include "./mesh2D.hh"

// functions for driving a simple post processing (printer,writer,postscript output)
#include "./mesh2D_printer.hh"
#include "./mesh2D_psout.hh"
#include "./mesh2D_rpars.hh"

const string ref_str[10] = { string("0"),  string("1"),  string("2"),  string("3"),
                             string("4"),  string("5"),  string("6"),  string("7"),
                             string("8"),  string("9") } ;

string change_blank( string inpstr ) {
  string retstr="" ;

  for ( int i=0; i<inpstr.length(); ++i ) {
    if ( inpstr[i]!=' ' ) {
      retstr += inpstr[i] ;
    } else {
      retstr += "-" ;
    }
  }
  return retstr ;
}

string get_filename( string inpstr ) {
  string retstr="" ;
  for ( int i=0; i<inpstr.length(); ++i ) {
    if ( inpstr[i]=='/' ) {
      retstr = "" ;
    } else {
      retstr += inpstr[i] ;
    }
  }
  return retstr ;
}

template<class MESH>
//void post_proc( int ilev, mesh_2Dv & mesh, Read_Mesh2D_Pars & rpar ) {
void post_proc( int ilev, MESH & mesh, Read_Mesh2D_Pars & rpar ) {

  string mesh_name  = get_filename(mesh.get_mesh_name()) ;
  string psout_file = string("./MESH-PS/") ;
  string outpt_file = string("./MESH-OUTP/") ;
  
  if ( ilev<10 ) { 
    mesh_name  += string("-Ref_")   + ref_str[ilev] ;
    psout_file += string("mesh2D_") + ref_str[ilev] ;
    outpt_file += string("mesh2D_") + ref_str[ilev] ;
  } else {
    mesh_name  += string("-Ref_")   + ref_str[ilev/10] + ref_str[ilev%10] ;
    psout_file += string("mesh2D_") + ref_str[ilev/10] + ref_str[ilev%10] ;
    outpt_file += string("mesh2D_") + ref_str[ilev/10] + ref_str[ilev%10] ;
  }
  
  // print contents of mesh database
  if ( ilev<=rpar.get_prtr_flag() ) {
    mesh2Dv_printer mesh_printer( mesh ) ;
    mesh_printer.print_all_datasets() ;
    mesh_printer.print_all_vertices() ;
    mesh_printer.print_all_faces   () ;
    mesh_printer.print_all_regions () ;
  }
  
  // print postscript picture of mesh
  if ( ilev<=rpar.get_psgd_flag() ) {

    //string psout_file = string("mesh") ;
    PS_Driver psd(psout_file) ;
    mesh2D_psout msh_psout( mesh, psd ) ;

    // just to change offset
    msh_psout.set_offset( OFFSET ) ;
    //msh_psout.set_offset( 1 ) ;
  
    //bool draw_bullets(false) ;
    bool draw_bullets = abs(rpar.get_numb_flag())/1000==1 ;
    bool draw_numbers = rpar.get_numb_flag()>=1 ;
    string mesh_name = mesh.get_mesh_name() ;
    if ( mesh_name==string("") ) { mesh_name=string("Mesh") ; }
    msh_psout . draw_mesh(mesh_name,draw_numbers) ;

    PRT( rpar.get_numb_flag() ) ;
    PRT( draw_bullets ) ;

    if ( draw_bullets ) {
      int flag = abs(rpar.get_numb_flag())-1000 ;
      PRT(flag) ;
      if ( flag/100      == 1 ) { msh_psout . draw_regn_bullet( 3 ) ; }
      if ( (flag/100)/10 == 1 ) { msh_psout . draw_face_bullet( 3 ) ; }
      if ( (flag/100)%10 == 1 ) { msh_psout . draw_vrtx_bullet( 3 ) ; }
    }
  }
}

template<class MESH>
void post_proc_meshless( int ilev, MESH & mesh, Read_Mesh2D_Pars & rpar ) {
  
  string mesh_name  = get_filename(mesh.get_mesh_name()) ;
  string psout_file = string("./MESH-PS/") ;
  string outpt_file = string("./MESH-OUTP/") ;
  
  if ( ilev<10 ) { 
    mesh_name  += string("-Ref_")   + ref_str[ilev] ;
    psout_file += string("mesh2D_") + ref_str[ilev] ;
    outpt_file += string("mesh2D_") + ref_str[ilev] ;
  } else {
    mesh_name  += string("-Ref_")   + ref_str[ilev/10] + ref_str[ilev%10] ;
    psout_file += string("mesh2D_") + ref_str[ilev/10] + ref_str[ilev%10] ;
    outpt_file += string("mesh2D_") + ref_str[ilev/10] + ref_str[ilev%10] ;
  }
  
  // print contents of mesh database
  if ( ilev<=rpar.get_prtr_flag() ) {
    mesh2Dv_printer mesh_printer( mesh ) ;
    mesh_printer.print_all_datasets() ;
    mesh_printer.print_all_vertices() ;
    mesh_printer.print_all_faces   () ;
    mesh_printer.print_all_regions () ;
  }
  
  // print postscript picture of mesh
  if ( ilev<=rpar.get_psgd_flag() ) {

    //string psout_file = string("mesh") ;
    PS_Driver psd(psout_file) ;
    mesh2D_psout msh_psout( mesh, psd ) ;

    // just to change offset
    msh_psout.set_offset( OFFSET ) ;
  
    // plot meshless points
    bool draw_bullets(true) ;
    bool draw_meshless(false) ;
    bool draw_numbers = rpar.get_numb_flag()>=1 ;
    string mesh_name = mesh.get_mesh_name() ;
    if ( mesh_name==string("") ) { mesh_name=string("Mesh") ; }
    msh_psout . draw_meshless(mesh_name,draw_numbers,draw_meshless) ;
  }
}

template<class MESH>
void post_proc_mesh( MESH & mesh, int n_msh ) {

  // print postscript picture of mesh
  string psout_file = string("./MESH-PS/mesh2D_") + ref_str[n_msh] ;
  string outpt_file = string("./MESH-OUTP/mesh2D") ;

  DBGF(BBB-00) ;
  PRT(psout_file) ;

  //string psout_file = string("mesh") ;
  PS_Driver psd(psout_file) ;
  mesh2D_psout msh_psout( mesh, psd ) ;

  // just to change offset
  msh_psout.set_offset( OFFSET ) ;
  //msh_psout.set_offset( 1 ) ;

  //bool draw_bullets(false) ;
  bool draw_bullets = false ;
  bool draw_numbers = false ;

  string mesh_name = mesh.get_mesh_name() ;
  if ( mesh_name==string("") ) { mesh_name=string("Mesh") ; }
  msh_psout . draw_mesh(mesh_name,draw_numbers) ;

  if ( draw_bullets ) {
    msh_psout . draw_regn_bullet( 3 ) ;
    msh_psout . draw_face_bullet( 3 ) ;
    msh_psout . draw_vrtx_bullet( 3 ) ;
  }

}

#endif // end of _MESH_2D_POST_HH
