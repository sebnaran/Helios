void mesh2Dv_builder::build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV,
				      vector<int> & regn_vrtx_list, vector<int> & regn_face_list,
				      vector<int> & face_vrtx_list, vector<int> & face_regn_list, 
				      vector<int> & fF ) {
  std::cout << "begin  ---> mesh2Dv_builder::build_the_mesh using all datasets <----" << std::endl << std::flush ;
  
  // check consistency
  assert( face_vrtx_list.size()==face_regn_list.size() ) ;
  assert( face_vrtx_list.size()%2==1 ) ;

  // -- set coordinates in mesh
  set_coords( xV, yV ) ;

  // determine nR 
  int nR = regn_vrtx_list.back() ;
  
  // -- core for building primary dataset
  build_RegnFace( regn_face_list, face_regn_list ) ;
  
  // set face-vertex data structure
  build_FaceVrtx( face_vrtx_list ) ;
  
  // -- set external flags for regions, faces and vertices
  mesh.fR.resize(nR) ;
  for ( int iR=0; iR<nR; ++iR ) { mesh.fR[iR] = UNSET ; }
  // ---
  int nF = face_vrtx_list.back() ;
  mesh.fF.resize(nF) ;
  for ( int iF=0; iF<nF; ++iF ) { mesh.fF[iF] = fF[iF] ; }
  // ---
  int nV = fV.size() ;
  mesh.fV.resize(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { mesh.fV[iV] = fV[iV] ; }  

  // -- tranpose datasets
  build_VrtxFace() ;
  build_FaceRegn() ;

  // -- build boundary lists
  build_boundary_lists() ;

  // -- compute/set last geometric quantities
  setup_geom_factors() ;
  
  // final reset of logical flags
  mesh.reset_status_flags() ;
  std::cout << "end of ---> mesh2Dv_builder::build_the_mesh using all datasets <----" << std::endl << std::flush ;
}

void mesh2Dv_builder :: build_RegnFace( vector<int> & regn_face_list,
                                        vector<int> & face_regn_list ) {
  std::cout << "begin build_RegnFace" << std::endl << std::flush ;
  int nR = regn_face_list.back() ;
  mesh.nR = nR ;
  mesh.RegnFace.setup( nR ) ;
  int k = 0 ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = regn_face_list[k++] ;
    mesh.RegnFace.setup(iR,nRF) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int  iF   = regn_face_list[k++] ;
      bool bval = iR==face_regn_list[2*iF] ? true : false ;
      mesh.RegnFace(iR,ilF) = iF ;
    }
  }
  std::cout << "end-->build_RegnFace" << std::endl << std::flush ;
}
