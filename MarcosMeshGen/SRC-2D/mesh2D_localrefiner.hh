#ifndef _MESH_2D_LOCALREFINER_HH
#define _MESH_2D_LOCALREFINER_HH

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

class LocalMeshRefinement {
private:
  struct sub_face {
    int face_id ;
    int flag_hz ;
    vector<int> sub_vrtx ;
  } ;

private:
  int nrx, nry ;
  int curr_iV, curr_iR ;
  mesh_2Dv & mesh ;
  vector<int> flag_cells, flag_faces ;
  vector<int> idx_bot, idx_top, idx_lft, idx_rgh ;
  vector<sub_face> split_face ;

  // new mesh input structure
  vector<int> regn_vlist, fR, fV ;
  vector<double> xV, yV ;

  bool is_face_vertical  ( int iF ) ;
  bool is_face_horizontal( int iF ) ;

  void build_face_substructure() ;
  void build_cell_substructure() ;
  void build_cell_substructure( int iR ) ;
  void build_cell_frame       ( int iR ) ;

  void build_unrefined_cells  () ;
  void build_unrefined_cells( int iR ) ;

public:
  LocalMeshRefinement( mesh_2Dv & _mesh ) : mesh(_mesh), nrx(0), nry(0) {} ;
  ~LocalMeshRefinement() {} ;

  void setup( vector<int> & flagged_cells_list, int _nrx, int _nry ) ;

  void build_local_mesh_refinement( mesh_2Dv & new_mesh ) {
    // ----------------------------------------
    build_face_substructure() ;
    build_cell_substructure() ;
    build_unrefined_cells  () ;
    //-----------------------------------------------------------------------------------------
    // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
    //-----------------------------------------------------------------------------------------
    mesh2Dv_builder mesh_builder(new_mesh) ;
    mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
    // set mesh name (standard flag 1xx)
    char str_nrx[10] ; sprintf(str_nrx, "%d",nrx) ;
    char str_nry[10] ; sprintf(str_nry, "%d",nry) ;
    string mesh_name = string("Local Mesh Refinement")+string("-")+str_nrx+string("x")+str_nry ;
    mesh . set_mesh_name( mesh_name ) ;
  }

} ;

bool LocalMeshRefinement::is_face_vertical( int iF ) {
  double xV0 = mesh.coords_V( mesh.face_vrtx(iF,0), 0 ) ;
  double xV1 = mesh.coords_V( mesh.face_vrtx(iF,1), 0 ) ;
  return abs(xV0-xV1)<1e-12 ;
}
bool LocalMeshRefinement::is_face_horizontal( int iF ) {
  double yV0 = mesh.coords_V( mesh.face_vrtx(iF,0), 1 ) ;
  double yV1 = mesh.coords_V( mesh.face_vrtx(iF,1), 1 ) ;
  return abs(yV0-yV1)<1e-12 ;
}

void LocalMeshRefinement::build_face_substructure() {
  MSGF("begin  build_face_substructure") ;

  // size
  int nV = mesh.n_vertex() ;
  int nF = mesh.n_face() ;

  // set curr_iV
  curr_iV = nV ;

  // set split_face
  split_face.resize(nF) ;
  for ( int iF=0; iF<nF; ++iF ) {
    if ( flag_faces[iF]==1 ) {
      split_face[iF].face_id = iF ;

      // VERTICAL FACE
      if ( is_face_vertical(iF) ) {

	split_face[iF].flag_hz = 0 ;
	split_face[iF].sub_vrtx.resize(nry+2) ;
	split_face[iF].sub_vrtx[0]     = mesh.face_vrtx(iF,0) ;
	split_face[iF].sub_vrtx[nry+1] = mesh.face_vrtx(iF,1) ;

	for ( int ilV=1; ilV<=nry; ++ilV ) {
	  split_face[iF].sub_vrtx[ilV] = curr_iV++ ;
	}

	double xface[2] ;
	double yface[2] ;
	for ( int il=0; il<2; ++il ) {
	  int iV = mesh.face_vrtx(iF,il) ;
	  xface[il] = mesh.coords_V(iV,0) ;
	  yface[il] = mesh.coords_V(iV,1) ;
	}
	double dx = (xface[1]-xface[0])/double(nry+1) ;
	double dy = (yface[1]-yface[0])/double(nry+1) ;
	for ( int ilV=1; ilV<=nry; ++ilV ) {
	  double xref = xface[0] + ilV*dx ;
	  double yref = yface[0] + ilV*dy ;
	  xV.push_back( xref  ) ;
	  yV.push_back( yref  ) ;
	  fV.push_back( UNSET ) ;
	}
	
	// HORIZONTAL FACE
      } else if ( is_face_horizontal(iF) ) {

	split_face[iF].flag_hz = 1 ;
	split_face[iF].sub_vrtx.resize(nrx+2) ;
	split_face[iF].sub_vrtx[0]     = mesh.face_vrtx(iF,0) ;
	split_face[iF].sub_vrtx[nrx+1] = mesh.face_vrtx(iF,1) ;

	for ( int ilV=1; ilV<=nrx; ++ilV ) {
	  split_face[iF].sub_vrtx[ilV] = curr_iV++ ;
	}

	double xface[2] ;
	double yface[2] ;
	for ( int il=0; il<2; ++il ) {
	  int iV = mesh.face_vrtx(iF,il) ;
	  xface[il] = mesh.coords_V(iV,0) ;
	  yface[il] = mesh.coords_V(iV,1) ;
	}
	double dx = (xface[1]-xface[0])/double(nrx+1) ;
	double dy = (yface[1]-yface[0])/double(nrx+1) ;
	for ( int ilV=1; ilV<=nrx; ++ilV ) {
	  double xref = xface[0] + ilV*dx ;
	  double yref = yface[0] + ilV*dy ;
	  xV.push_back( xref  ) ;
	  yV.push_back( yref  ) ;
	  fV.push_back( UNSET ) ;
	}

      } else {
	MSGF("the mesh is not orthogonal!") ;
	assert(false) ;
      }
    }
  } // end of -->> for ( int iF=0; iF<nF; ++iF ) {...
    MSGF("end of build_face_substructure") ;
}

void LocalMeshRefinement::build_unrefined_cells() {
  MSGF("begin  build_unrefined_cells") ;
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    if ( flag_cells[iR]==0 ) {
      build_unrefined_cells( iR ) ;
    }
  }
  // set final number of regions in the refined mesh
  regn_vlist.push_back( curr_iR ) ;
  MSGF("end of build_cell_substructure") ;
}

void LocalMeshRefinement::build_unrefined_cells( int iR ) {

  // set a temporary list of vertices
  vector<int> cell_vrtx ;

  int nRF = mesh.n_regn_face(iR) ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = mesh.regn_face(iR,ilF) ;
    if ( flag_faces[iF]==0 ) {
      if ( mesh.ok_regn_face(iR,ilF) ) {
	cell_vrtx.push_back( mesh.face_vrtx(iF,0) ) ;
      } else { 
	cell_vrtx.push_back( mesh.face_vrtx(iF,1) ) ;
      }
    } else {
      int nsize= split_face[iF].sub_vrtx.size() ;
      if ( mesh.ok_regn_face(iR,ilF) ) {
	for ( int ilV=0; ilV<nsize-1; ++ilV ) {
	  cell_vrtx.push_back( split_face[iF].sub_vrtx[ilV] ) ;
	}
      } else {
	for ( int ilV=0; ilV<nsize-1; ++ilV ) {
	  cell_vrtx.push_back( split_face[iF].sub_vrtx[nsize-1-ilV] ) ;
	} 
      }
    } // end of -->> if ( flag_faces[iF]==0 ) {...
  } // end of -->> for ( int ilF=0; ilF<nRF; ++ilF ) {...

  // set a new region
  curr_iR++ ;
  regn_vlist.push_back( cell_vrtx.size() ) ;
  for ( int il=0; il<cell_vrtx.size(); ++il ) {
    regn_vlist.push_back( cell_vrtx[il] ) ;
  }
  fR.push_back( UNSET ) ;
    
}

void LocalMeshRefinement::build_cell_substructure() {
  MSGF("begin  build_cell_substructure") ;
  int nR = mesh.n_region() ;
  curr_iR = 0 ;
  for ( int iR=0; iR<nR; ++iR ) {
    if ( flag_cells[iR]==1 ) {
      build_cell_substructure( iR ) ;
    }
  }
  MSGF("end of build_cell_substructure") ;
}

void LocalMeshRefinement::build_cell_substructure( int iR ) {
  //MSGF("begin  build_cell_substructure: iR = "<<iR) ;

  // build frame data structure
  build_cell_frame( iR ) ;

  // build full cell substructure
  int cell_vrtx[nrx+2][nry+2] ;

  // set UNSET status
  for ( int ilr=0; ilr<nrx+2; ++ilr ) {
    for ( int jlr=0; jlr<nry+2; ++jlr ) {
      cell_vrtx[ilr][jlr] = UNSET ;
    }
  }

  // set left & right
  for ( int jlr=0; jlr<nry+2; ++jlr ) {
    cell_vrtx[0]    [jlr] = idx_lft[jlr] ;
    cell_vrtx[nrx+1][jlr] = idx_rgh[jlr] ;
  }

  // set top & bottom
  for ( int ilr=0; ilr<nrx+2; ++ilr ) {
    cell_vrtx[ilr][0]     = idx_bot[ilr] ;
    cell_vrtx[ilr][nry+1] = idx_top[ilr] ;
  }

  // add internal vertices
  double x0 = mesh.coords_V( idx_bot[0], 0 ) ;
  double y0 = mesh.coords_V( idx_bot[0], 1 ) ;
  double dx = ( mesh.coords_V(idx_rgh[0],0)-mesh.coords_V(idx_lft[0],0) )/double(nrx+1) ; 
  double dy = ( mesh.coords_V(idx_top[0],1)-mesh.coords_V(idx_bot[0],1) )/double(nry+1) ; 
  for ( int ilr=0; ilr<nrx+2; ++ilr ) {
    for ( int jlr=0; jlr<nry+2; ++jlr ) {
      if ( cell_vrtx[ilr][jlr] == UNSET ) {
	cell_vrtx[ilr][jlr] = curr_iV++ ;
	double xref = x0 + ilr*dx ;
	double yref = y0 + jlr*dy ;
	xV.push_back( xref  ) ;
	yV.push_back( yref  ) ;
	fV.push_back( UNSET ) ;
      } // END OF -->> if ( cell_vrtx[ilr][jlr] == UNSET ) {...
    }
  }

  // build the internal substructure
  for ( int jlr=0; jlr<nry+1; ++jlr ) {    
    for ( int ilr=0; ilr<nrx+1; ++ilr ) {
      // set a new region: (i,j), (i+1,j), (i+1,j+1), (i,j+1)
      curr_iR++ ;
      regn_vlist.push_back( 4 ) ;
      regn_vlist.push_back( cell_vrtx[ilr  ][jlr  ] ) ;
      regn_vlist.push_back( cell_vrtx[ilr+1][jlr  ] ) ;
      regn_vlist.push_back( cell_vrtx[ilr+1][jlr+1] ) ;
      regn_vlist.push_back( cell_vrtx[ilr  ][jlr+1] ) ;
      fR.push_back( UNSET ) ;
    }
  }
  //MSGF("end of build_cell_substructure: iR = "<<iR) ;
}

void LocalMeshRefinement::build_cell_frame( int iR ) {
  //MSGF("begin  build_cell_frame: iR = "<<iR) ;

  idx_bot.resize(0) ;
  idx_top.resize(0) ;
  idx_lft.resize(0) ;
  idx_rgh.resize(0) ;

  int nRF = mesh.n_regn_face(iR) ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {

    int iF = mesh.regn_face(iR,ilF) ;
    int nr = split_face[iF].sub_vrtx.size() ;
    int hz = split_face[iF].flag_hz ;

    //VAL(ilF) ; PRT(iF) ;

    // horizontal faces
    if ( hz==1 ) {

      if ( idx_bot.size()==0 ) {
	for ( int ilr=0; ilr<nr; ++ilr ) {
	  if ( mesh.ok_regn_face(iR,ilF) ) {
	    idx_bot.push_back( split_face[iF].sub_vrtx[ilr] ) ;
	  } else {
	    idx_bot.push_back( split_face[iF].sub_vrtx[nr-1-ilr] ) ;
	  }
	}
      } else {
	if ( mesh.ok_regn_face(iR,ilF) ) {
	  for ( int ilr=0; ilr<nr; ++ilr ) {
	    idx_top.push_back( split_face[iF].sub_vrtx[nr-1-ilr] ) ;
	  }
	} else {
	  for ( int ilr=0; ilr<nr; ++ilr ) {
	    idx_top.push_back( split_face[iF].sub_vrtx[ilr] ) ;
	  }
	}
      }

    } // end of -->> if ( hz==1 ) {...

    // vertical faces
    if ( hz==0 ) {
      if ( idx_rgh.size()==0 ) {
	for ( int ilr=0; ilr<nr; ++ilr ) {
	  if ( mesh.ok_regn_face(iR,ilF) ) {
	    idx_rgh.push_back( split_face[iF].sub_vrtx[ilr] ) ;
	  } else {
	    idx_rgh.push_back( split_face[iF].sub_vrtx[nr-1-ilr] ) ;
	  }
	}
      } else {
	if ( mesh.ok_regn_face(iR,ilF) ) {
	  for ( int ilr=0; ilr<nr; ++ilr ) {
	    idx_lft.push_back( split_face[iF].sub_vrtx[nr-1-ilr] ) ;
	  }
	} else {
	  for ( int ilr=0; ilr<nr; ++ilr ) {
	    idx_lft.push_back( split_face[iF].sub_vrtx[ilr] ) ;
	  }
	}
      }
    } // end of -->> if ( hz==o ) {...
  } // end of for ( int ilF=0; ilF<nRF; ++ilF ) {...

  // swap idx_bot and idx_top if we need it
  if ( mesh.coords_V(idx_bot[0],1)>mesh.coords_V(idx_top[0],1) ) {
    swap( idx_bot, idx_top ) ;
  }

  // swap idx_bot and idx_top if we need it
  if ( mesh.coords_V(idx_lft[0],0)>mesh.coords_V(idx_rgh[0],0) ) {
    swap( idx_lft, idx_rgh ) ;
  }

  //MSGF("end of build_cell_frame: iR = "<<iR) ;
}

void LocalMeshRefinement::setup( vector<int> & flagged_cells_list, int _nrx, int _nry ) {

  // set input refinement size
  nrx = _nrx ;
  nry = _nry ;

  // set flag_cells, flag_faces
  int nR = mesh.n_region() ;
  int nV = mesh.n_vertex() ;
  int nF = mesh.n_face  () ;

  // resize flag_cells and flag_faces arrays
  flag_cells.resize(nR) ;
  flag_faces.resize(nF) ;

  // set default value --->>> 0
  for ( int iR=0; iR<nR; ++iR ) { flag_cells[iR] = 0 ; }
  for ( int iF=0; iF<nF; ++iF ) { flag_faces[iF] = 0 ; }
    
  // set flagged cells --->>> 1
  int nfR = flagged_cells_list.size() ;
  for ( int ilR=0; ilR<nfR; ++ilR ) {
    int iR = flagged_cells_list[ilR] ;
    flag_cells[iR] = 1 ;
  }

  // set flagged faces
  for ( int iF=0; iF<nF; ++iF ) {
    int iR0 = mesh.face_regn( iF, 0 ) ;
    int iR1 = mesh.face_regn( iF, 1 ) ;
    if ( mesh.is_internal_face(iF) ) {
      if ( flag_cells[iR0]==1 || flag_cells[iR1]==1 ) { flag_faces[iF] = 1 ; } 
    } else {
      if ( flag_cells[iR0]==1 ) { flag_faces[iF] = 1 ; } 
    }
  }

#if 0
  // TBR
  PRT_ARR( flag_cells ) ;
  LINE(--) ;
  PRT_ARR( flag_faces ) ;
  LINE(--) ;
  exit(0) ;
#endif

  // initialize the mesh
  xV.resize(nV) ;
  yV.resize(nV) ;
  fV.resize(nV) ;
  for ( int iV=0; iV<nV; ++iV ) {
    xV[iV] = mesh.coords_V(iV,0) ;
    yV[iV] = mesh.coords_V(iV,1) ;
    fV[iV] = UNSET ;
  }

}

#endif // end of _MESH_2D_LOCALREFINER_HH

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

#if 0

// EXAMPLE

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
  flagged_cells_list.push_back(0) ;
  flagged_cells_list.push_back(3) ;

  LocalMeshRefinement loc_msh_ref(mesh) ;
  loc_msh_ref.setup( flagged_cells_list ) ;

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

#endif

