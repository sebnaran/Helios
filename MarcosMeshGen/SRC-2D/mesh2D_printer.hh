class mesh2Dv_printer {

  static const int DIM = 2 ;

private:
  mesh_2Dv & mesh ;
  int offset ;

public:
  mesh2Dv_printer( mesh_2Dv & _mesh, int _offset=0 ) : mesh(_mesh), offset(_offset) {}
  ~mesh2Dv_printer() {}

  // print DATASETS
  void print_FaceVrtx( ostream & LOGF ) ;
  void print_RegnFace( ostream & LOGF ) ;
  // ---
  void print_FaceRegn( ostream & LOGF ) ;
  void print_VrtxFace( ostream & LOGF ) ;
  // ---
  void print_boundary_lists( ostream & LOGF ) ;
  // --- 
  void print_all_datasets() ;

  // print REGIONS
  void print_all_regions() ;
  void print_region( int iR, ostream & LOGF ) ;

  // print FACES
  void print_all_faces() ;
  void print_face( int iF, ostream & LOGF ) ;

  // print VERTICES
  void print_all_vertices() ;
  void print_vertex( int iV, ostream & LOGF ) ;

  // print MESH additional data
  void print_additional_data() ;
} ;

// --------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_FaceVrtx( ostream & LOGF ) {
  // FaceVrtx
  const int nFV = 2 ;
  LOGF << "FaceVrtx " << endl ;
  LOGF << "#of faces = "  << mesh.n_face() << endl ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    LOGF << "iF = " << iF+offset << ", { " 
	 << mesh.FaceVrtx[nFV*iF+0]+offset << ", " 
	 << mesh.FaceVrtx[nFV*iF+1]+offset
	 << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh2Dv_printer :: print_RegnFace( ostream & LOGF ) {
  // RegnFace
  LOGF << "RegnFace " << endl ;
  LOGF << "#of regns = "  << mesh.n_region() << endl ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    LOGF << "iR = " << iR+offset << ", { "  ;
    if ( mesh.RegnFace.size_loc(iR)>0 ) {
      for ( int j=0; j<mesh.RegnFace.size_loc(iR)-1; ++j ) {
	LOGF << "(" << mesh.ok_regn_face(iR,j) << ")" << mesh.RegnFace(iR,j)+offset << ", " ;
      }
      int j = mesh.RegnFace.size_loc(iR)-1 ;
	  LOGF << "(" << mesh.ok_regn_face(iR,j) << ")" 
	   << mesh.RegnFace(iR,j)+offset ;
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
// ------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_VrtxFace( ostream & LOGF ) {
  // VrtxFace
  LOGF << "VrtxFace " << endl ;
  LOGF << "#of vertices = "  << mesh.n_vertex() << flush << endl ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    LOGF << "iV = " << iV+offset << ", { " ;
    if ( mesh.VrtxFace.size_loc(iV)>0 ) {
      for ( int j=0; j<mesh.VrtxFace.size_loc(iV)-1; ++j ) {
	LOGF << "(" << mesh.ok_vrtx_face(iV,j) << ")" << mesh.VrtxFace(iV,j)+offset << ", " ;
      }
      int j = mesh.VrtxFace.size_loc(iV)-1 ;
      LOGF << "(" << mesh.ok_vrtx_face(iV,j) << ")" 
	   << mesh.VrtxFace(iV,j)+offset ; 
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh2Dv_printer :: print_FaceRegn( ostream & LOGF ) {
  // FaceRegn
  LOGF << "FaceRegn " << endl ;
  LOGF << "#of faces = "  << mesh.n_face() << endl ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    LOGF << "iF = " << iF+offset << ", { " 
	 << mesh.FaceRegn[2*iF+0]+offset << ", " << mesh.FaceRegn[2*iF+1]+offset 
	 << " }; flag = " << mesh.get_fF(iF) << endl ;
  }
  LOGF << "-------------------" << endl ;
}

void mesh2Dv_printer :: print_boundary_lists( ostream & LOGF ) {
  // boundary lists
  LOGF << "boundary lists " << endl ;
  LOGF << "#of boundary vertices = "  << mesh.n_bvertex() << endl ;
  for ( int i=0; i<mesh.n_bvertex(); ++i ) {
    LOGF << "bnd_vrtx[" << i << "] = " << mesh.bnd_vrtx[i]+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
  LOGF << "#of boundary faces = "  << mesh.n_bface() << endl ;
  for ( int i=0; i<mesh.n_bface(); ++i ) {
    LOGF << "bnd_face[" << i << "] = " << mesh.bnd_face[i]+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
  LOGF << "#of boundary regions = "  << mesh.n_bregion() << endl ;
  for ( int i=0; i<mesh.n_bregion(); ++i ) {
    LOGF << "bnd_regn[" << i << "] = " << mesh.bnd_regn[i]+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_all_datasets() {
  ofstream LOGF("dset.log") ;
  print_FaceVrtx( LOGF ) ;
  print_RegnFace( LOGF ) ;
  print_VrtxFace( LOGF ) ;
  print_FaceRegn( LOGF ) ;
  LOGF << "-------------------" << endl ;
  print_boundary_lists( LOGF ) ;
  LOGF.close() ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
# define print_list_offset(VEC)                                     \
  do {                                                       \
    LOGF << " { " ;					     \
    if ( VEC.size()>0 ) {				     \
      for ( int i=0 ; i<VEC.size()-1 ; ++i ) {		     \
	LOGF << VEC[i]+offset << ", " ;			     \
      }							     \
      LOGF << VEC.back()+offset ;				     \
    }							     \
    LOGF << " }" << endl ;				     \
  } while(0)
// --------------------------------------------------------------------------------------------
// REGIONS
// --------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_all_regions() {
  ofstream LOGF("region.log") ;
  LOGF << "#of regions: nR=" << mesh.n_region() << endl ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    print_region( iR, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh2Dv_printer :: print_region( int iR, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "region: iR =" << iR+offset << endl ;
  // faces
  int nRF = mesh.RegnFace.size_loc(iR);
  LOGF << "#faces: nRF=" << nRF << ", { " ;
  for ( int ilF=0; ilF<nRF-1; ++ilF ) {
    LOGF << "(" 
	 << mesh.ok_regn_face(iR,ilF) << ")" 
	 << mesh.regn_face   (iR,ilF)+offset << ", " ;
  }
  LOGF << "(" 
       << mesh.ok_regn_face(iR,nRF-1) << ")" 
       << mesh.regn_face   (iR,nRF-1)+offset << " }" << endl ;
  LOGF << flush ;
  // topological info
  vector<int> rlist, flist, vlist ;
  mesh.get_regn_regn( iR, rlist ) ;
  mesh.get_regn_face( iR, flist ) ;
  mesh.get_regn_vrtx( iR, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nRR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nRF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#vertices: nRV=" << vlist.size() ; print_list_offset(vlist) ;
  // geometrical info
  LOGF << "---" << endl ;
  LOGF << "#center: ( " 
       << mesh.coords_R(iR,0) << ", " << mesh.coords_R(iR,1) << " )" << endl << flush ; 
  LOGF << "#area:     " << mesh.get_regn_measure(iR) << endl << flush ; 
}
// --------------------------------------------------------------------------------------------
// FACES
// --------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_all_faces() {
  ofstream LOGF("face.log") ;
  LOGF << "#of faces: nF=" << mesh.n_face() << endl ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    print_face( iF, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh2Dv_printer :: print_face( int iF, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "face: iF =" << iF+offset << endl ;
  // topological info
  vector<int> rlist, flist, vlist ;
  mesh.get_face_regn( iF, rlist ) ;
  mesh.get_face_face( iF, flist ) ;
  mesh.get_face_vrtx( iF, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nFR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nFF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#vertices: nFV=" << vlist.size() ; print_list_offset(vlist) ;
  // geometrical info
  LOGF << "---" << endl ;
  LOGF << "#center: ( " 
       << mesh.coords_F(iF,0) << ", " 
       << mesh.coords_F(iF,1) << " )" << endl ; 
  LOGF << "#nrm: ( " 
       << mesh.get_nor(iF,0) << ", " 
       << mesh.get_nor(iF,1) << " )" << endl ; 
  LOGF << "#tng: ( " 
       << mesh.get_tng(iF,0) << ", " 
       << mesh.get_tng(iF,1) << " )" << endl ; 
  LOGF << "#length:   " << mesh.get_face_measure(iF) << endl ; 
}
// --------------------------------------------------------------------------------------------
// VERTICES
// --------------------------------------------------------------------------------------------
void mesh2Dv_printer :: print_all_vertices() {
  ofstream LOGF("vertex.log") ;
  LOGF << "#of vertices: nV=" << mesh.n_vertex() << endl ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    print_vertex( iV, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh2Dv_printer :: print_vertex( int iV, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "vertex: iV =" << iV+offset << endl ;
  // coordinates
  LOGF << "#coords:= ( " 
       << mesh.vrtx_coords[DIM*iV+0] << ", " 
       << mesh.vrtx_coords[DIM*iV+1] << ") " << endl << flush ; 
  // topological info
  vector<int> rlist, flist, vlist ;
  mesh.get_vrtx_regn( iV, rlist ) ;
  mesh.get_vrtx_face( iV, flist ) ;
  mesh.get_vrtx_vrtx( iV, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nVR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nVF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#vertices: nVV=" << vlist.size() ; print_list_offset(vlist) ;
}
void mesh2Dv_printer :: print_additional_data() {
  ofstream LOGF("mesh.log") ;
  double xmin, ymin, xmax, ymax ;  
  mesh.bbox( xmin, ymin, xmax, ymax  ) ;
  LOGF << "xmin = " << xmin << " mesh.xmin() = " << mesh.xmin() << endl ;
  LOGF << "ymin = " << ymin << " mesh.ymin() = " << mesh.ymin() << endl ;
  LOGF << "xmax = " << xmax << " mesh.xmax() = " << mesh.xmax() << endl ;
  LOGF << "ymax = " << ymax << " mesh.ymax() = " << mesh.ymax() << endl ;
  LOGF << "mesh.h_max() = " << mesh.h_max() ;
  LOGF << "mesh.h_min() = " << mesh.h_min() ;
  LOGF << "mesh.h_sqr() = " << mesh.h_sqr() ;
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    LOGF << "iR = " << iR 
	 << "mesh.h_regn(iR) = "      << mesh.h_regn(iR) 
	 << "mesh.h_regn_face(iR) = " << mesh.h_regn_face(iR) 
	 << endl ;
  }
  LOGF.close() ;
} 
