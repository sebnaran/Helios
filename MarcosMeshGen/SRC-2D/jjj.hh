//-----------------------------------------------------------------------------------------
void build_dual_mesh( mesh_2Dv & mesh, mesh_2Dv & dual_mesh ) {
  DBGF(START build_dual_mesh) ;
  int nV = mesh.n_region()+mesh.n_bface()+mesh.n_bvertex() ; // #vertices of dual mesh
  int nR = mesh.n_vertex() ;                                 // #regions  of dual mesh
  assert( nV>0 && nR>0 ) ;
  vector<double> xV, yV ;
  vector<int>    regn_vlist ;
  vector<int>    fV(nV), fR(nR) ;
  for ( int iV=0; iV<nV; ++iV ) { fV[iV]=UNSET ; }
  for ( int iR=0; iR<nR; ++iR ) { fR[iR]=UNSET ; }
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF DUAL MESH
  //-----------------------------------------------------------------------------------------
  // dual vrtx == region center
  vector<int> face_idx(mesh.n_face(),  UNSET) ;
  vector<int> vrtx_idx(mesh.n_vertex(),UNSET) ;
  for ( int iR=0 ; iR<mesh.n_region() ; ++iR ) {
    xV.push_back( mesh.coords_R(iR,0) ) ;
    yV.push_back( mesh.coords_R(iR,1) ) ;
  }
  int new_iV = mesh.n_region() ;
  // dual vrtx on bnd == face midpoint 
  for ( int ilF=0 ; ilF<mesh.n_bface() ; ++ilF ) {
    int ibF = mesh.get_bnd_face(ilF) ;
    xV.push_back( mesh.coords_F(ibF,0) ) ;
    yV.push_back( mesh.coords_F(ibF,1) ) ;
    face_idx[ibF] = new_iV++ ;
  }
  // dual vrtx on bnd == bnd vrtx
  for ( int ilV=0 ; ilV<mesh.n_bvertex() ; ++ilV ) {
    int ibV = mesh.get_bnd_vrtx(ilV) ;
    xV.push_back( mesh.coords_V(ibV,0) ) ;
    yV.push_back( mesh.coords_V(ibV,1) ) ;
    vrtx_idx[ibV] = new_iV++ ;
  }
  assert( nV==xV.size() && xV.size()==yV.size() && new_iV==nV ) ;
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF DUAL MESH
  // set region vertex lists
  // questa routine e` molto piu` complicata dell'altra perche` non assume
  // un ordinamento forte del dataset vrtx_face
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    vector<int> flist ;
    mesh.get_vrtx_face( iV, flist ) ;
    int nVF = flist.size() ;

    vector<bool> mask_F(nVF,true) ;

    int next_iR = UNSET ;
    int iF      = UNSET ;
    int iR0(UNSET), iR1(UNSET) ;

    if ( mesh.is_boundary_vrtx(iV) ) {
      
      bool searching_first_face = true ;
      for ( int ilF=0; ilF<nVF && searching_first_face; ++ilF ) {
	iF = flist[ilF] ;
	if ( mesh.face_regn(iF,1)==UNSET && mesh.face_vrtx(iF,0)==iV ) {
	  searching_first_face = false ;
	  mask_F[ ilF ] = false; 
	} // if (true) ==> out with current value of iF
      }

      int iR0 = mesh.face_regn( iF, 0 ) ;
      
      regn_vlist.push_back( nVF+2 ) ;
      regn_vlist.push_back( vrtx_idx[iV] ) ;
      regn_vlist.push_back( face_idx[iF] ) ;
      regn_vlist.push_back( iR0 ) ;
      
      next_iR = iR0 ;

      assert( vrtx_idx[iV]!=UNSET ) ;
      assert( face_idx[iF]!=UNSET ) ;
      assert( mesh.is_boundary_face(iF) ) ;

    } else {
      
      // start sequence
      regn_vlist.push_back( nVF ) ;

      iF  = flist[ 0 ] ;
      mask_F[ 0 ]  = false; 
      
      iR0 = mesh.face_regn( iF, 0 ) ;
      iR1 = mesh.face_regn( iF, 1 ) ;

      if ( iV==mesh.face_vrtx(iF,0) ) { // F is coming out of V
	// regn sequence: iR1, iR0 (counterclock-wise orientation)
	//regn_vlist.push_back( iR1 ) ;
	regn_vlist.push_back( iR0 ) ;
	next_iR = iR0 ;
      } else {
	// regn sequence: iR0, iR1 (counterclock-wise orientation)
	//regn_vlist.push_back( iR0 ) ;
	regn_vlist.push_back( iR1 ) ;
	next_iR = iR1 ;
      }
    }
      
    // iF is set from previous loop
    for ( int ilF=0; ilF<nVF; ++ilF ) {
      bool ok_found = search_for_next_regn( iF, next_iR, mesh, flist, mask_F ) ;
      if ( ok_found ) {
	if ( next_iR!=UNSET ) { 
	  regn_vlist.push_back( next_iR ) ;
	} else {
	  assert( mesh.is_boundary_face(iF) ) ;
	  assert( ilF==nVF-2 ) ;
	}
      }
    }
    for ( int ilF=0; ilF<nVF; ++ilF ) { // check everything is OK
      assert( !mask_F[ilF] ) ;
    }

    if ( mesh.is_boundary_vrtx(iV) ) { // if bnd vertex, set last face mid-point
      assert( mesh.is_boundary_face(iF) ) ;
      assert( face_idx[iF]!=UNSET ) ;
      regn_vlist.push_back( face_idx[iF] ) ;
    }
  }
  regn_vlist.push_back( nR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR DUAL MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder dual_mesh_builder(dual_mesh) ;
  dual_mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;

  string mesh_name = string("Dual-of-") + mesh.get_mesh_name() ;
  dual_mesh.set_mesh_name( mesh_name ) ;
  DBGF(END OF build_dual_mesh) ;
}
