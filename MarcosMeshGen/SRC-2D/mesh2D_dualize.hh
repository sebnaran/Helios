#ifndef _MESH_2D_DUALIZE_HH
#define _MESH_2D_DUALIZE_HH

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// DUAL MESH CONSTRUCTION:
// - it builds a dual mesh from a domain partition of CONVEX polygons;
// - if the starting mesh is a Delaunay triangulation, the dual mesh is (quasi-)Voronoi;
// - external flags are set to UNSET
// - if boolean true_voronoi==true then it uses circumcenters, otherwise barycenters 
// - using: mesh_2Dv<..> and mesh2D_builder<..> 
// 
// NOTE: if the starting mesh is non composed by CONVEX cells, the result is unpredictable
//-----------------------------------------------------------------------------------------

// input:   next_iR (current), mesh, flist, mask_F
// output:  next_iF, next_iR [new if available, otherwise UNSET] 
bool search_for_next_regn( int          & next_iF, 
			   int          & next_iR,
			   mesh_2Dv     & mesh, 
			   vector<int>  & flist, 
			   vector<bool> & mask_F ) {
  bool ok_found = false ;
  for ( int ilF=0; ilF<flist.size() && !ok_found; ++ilF ) {
    if ( mask_F[ilF] ) {
      int iF  = flist[ilF] ;
      int iR0 = mesh.face_regn(iF,0) ;
      int iR1 = mesh.face_regn(iF,1) ;
      if ( next_iR == iR0 || next_iR == iR1 ) { // !!! found a new face !!!
	next_iF = iF ;
	next_iR = next_iR==iR0 ? iR1 : iR0 ;
	mask_F[ilF] = false ;
	ok_found    = true  ;
      }
    }
  }
  return ok_found ;
}

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

#include "pqq.hh"

bool vertex_is_an_hanging_node_of_dual_mesh( int iV, mesh_2Dv & mesh ) {
  bool retval = false ;
  if ( mesh.is_boundary_vrtx(iV) ) {
    
    // get the list of the faces incident iV
    vector<int> flist ;
    mesh.get_vrtx_face( iV, flist ) ;
    int nVF = flist.size() ;
    
    // get the two boundary faces
    int ibF0(UNSET), ibF1(UNSET) ;
    for ( int ilF=0; ilF<nVF; ++ilF ) {
      int iF = flist[ilF] ;
      if ( mesh.is_boundary_face(iF) ) {
	if ( ibF0==UNSET ) { 
	  ibF0=iF ; 
	} else {
	  ibF1=iF ;
	  break ;
	}
      }
    }

    // check we got the faces
    assert( ibF0!=UNSET && ibF1!=UNSET ) ;

    // check if the normal are parallel
    double nxF0 = mesh.get_nor( ibF0, 0 ) ;
    double nyF0 = mesh.get_nor( ibF0, 1 ) ;
    double nxF1 = mesh.get_nor( ibF1, 0 ) ;
    double nyF1 = mesh.get_nor( ibF1, 1 ) ;
    if ( abs( 1.-nxF0*nxF1-nyF0*nyF1 )<1.e-8 ) {
      retval = true ;
    }
  }
  return retval ;
}

//-----------------------------------------------------------------------------------------
// determine if the boundary vertex is a potential hanging node in the dual mesh
// note: internal vertices are not processed and the function returns false because a 
// vertex can never be a hanging node
void build_dual_mesh_no_boundary_hanging_nodes( mesh_2Dv & mesh, mesh_2Dv & dual_mesh ) {
  DBGF(START build_dual_mesh) ;
  vector<double> xV, yV ;
  vector<int>    regn_vlist, fV, fR ;
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
  //dual vrtx on bnd == face midpoint
  for ( int ilF=0 ; ilF<mesh.n_bface() ; ++ilF ) {
    int ibF = mesh.get_bnd_face(ilF) ;
    xV.push_back( mesh.coords_F(ibF,0) ) ;
    yV.push_back( mesh.coords_F(ibF,1) ) ;
    face_idx[ibF] = new_iV++ ;
  }
  // dual vrtx on bnd == bnd vrtx
  for ( int ilV=0 ; ilV<mesh.n_bvertex() ; ++ilV ) {
    int ibV = mesh.get_bnd_vrtx(ilV) ;
    if ( !vertex_is_an_hanging_node_of_dual_mesh( ibV, mesh ) ) {
      xV.push_back( mesh.coords_V(ibV,0) ) ;
      yV.push_back( mesh.coords_V(ibV,1) ) ;
      vrtx_idx[ibV] = new_iV++ ;
    }
  }
  assert( new_iV==xV.size() && xV.size()==yV.size() ) ;
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF DUAL MESH
  // set region vertex lists
  // questa routine e` molto piu` complicata dell'altra perche` non assume
  // un ordinamento forte del dataset vrtx_face
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {

    // build region's flag
    fR.push_back( UNSET ) ;

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

      if ( vrtx_idx[iV]!=UNSET ) {
	regn_vlist.push_back( nVF+2 ) ;
	regn_vlist.push_back( vrtx_idx[iV] ) ;
      } else {
	regn_vlist.push_back( nVF+1 ) ;
      }
      regn_vlist.push_back( face_idx[iF] ) ;
      regn_vlist.push_back( iR0 ) ;
      
      next_iR = iR0 ;

      //assert( vrtx_idx[iV]!=UNSET ) ;
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
  int nR = mesh.n_vertex() ;
  assert( fR.size()==nR ) ;
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

//-----------------------------------------------------------------------------------------
// this routine gets rid of the hanging nodes on the boundary edges
// input:  mesh (already built) with possible hanging nodes
// output: mesh (use the same reference) without hanging nodes
//-----------------------------------------------------------------------------------------
void clean_up_boundary_nodes( mesh_2Dv & mesh ) {
  MSGF("begin  clean_up_boundary_nodes") ;

  vector<double> xV, yV ;
  vector<int> fV, fR, regn_vlist ;

  int new_iV = 0 ;
  int nV = mesh.n_vertex() ;
  vector<int> new_vrtx_idx(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { 

    // default 
    bool ok_node = true ;

    // check if the node must be removed from the mesh
    if ( mesh.is_boundary_vrtx(iV) ) {
      vector<int> vrtx_flist ;      
      mesh.get_vrtx_face( iV, vrtx_flist ) ;
      //VAL(iV) ; PRT( vrtx_flist.size() ) ;
      if ( vrtx_flist.size()==2 ) {
	int iF0 = vrtx_flist[0] ;
	int iF1 = vrtx_flist[1] ;
	double dnx = abs(mesh.get_nor(iF1,0)-mesh.get_nor(iF0,0)) ;
	double dny = abs(mesh.get_nor(iF1,1)-mesh.get_nor(iF0,1)) ;
	if ( dnx+dny<1e-12 ) {
	  ok_node = false ;
	  //cout << "<<<--------" ;
	}
	//cout << endl ;
      }
    } // end of --->>> if ( mesh.is_boundary_vrtx(iV) ) {...
    
    if ( ok_node ) {
      new_vrtx_idx[iV] = new_iV++ ;
      xV.push_back( mesh.coords_V(iV,0) ) ;
      yV.push_back( mesh.coords_V(iV,1) ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
    } else {
      new_vrtx_idx[iV] = UNSET ;
    }
  }
  
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA FOR SPLIT MESH (final construction)
  //-----------------------------------------------------------------------------------------

  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    
    // get vertex list
    vector<int> rvlist ;
    mesh.get_regn_vrtx( iR, rvlist ) ;
    int nRV = rvlist.size() ;

    int new_nRV = 0 ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = rvlist[ilV] ;
      if ( new_vrtx_idx[iV] != UNSET ) {
	new_nRV++ ;
      }
    }
   
    regn_vlist.push_back( new_nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = rvlist[ilV] ;
      if ( new_vrtx_idx[iV] != UNSET ) {
	regn_vlist.push_back( new_vrtx_idx[iV] ) ;
      }
    }
    fR.push_back( mesh.get_fR(iR) ) ;
  }

  regn_vlist.push_back( nR ) ;

  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR SPLIT MESH (final construction)
  //-----------------------------------------------------------------------------------------
  // build new mesh
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;

  MSGF("end of clean_up_boundary_nodes") ;
}


#endif // end of _MESH_2D_DUALIZE_HH
