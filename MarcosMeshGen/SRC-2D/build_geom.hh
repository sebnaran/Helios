// build geometric quantities
void mesh2Dv_builder :: setup_geom_factors() {
  // setup of mesh arrays
  // regions:
  int nR = mesh.n_region() ;
  mesh.regn_area.resize(nR) ;
  mesh.regn_coords.resize(nR*DIM) ;
  // faces:
  int nF = mesh.n_face() ;
  mesh.face_length.resize(nF) ;
  // ----------------------
#if 1
  set_regn_geom_factors() ;
  set_face_geom_factors() ;
#else
    cerr << ">>>>>>>>>>>REMOVED GEOMETRIC FACTORS!!!<<<<<<<<<<<<<<" << endl << flush ;
#endif
}
void mesh2Dv_builder :: set_regn_geom_factors() { 
  // this routine works nicely if R is star-shaped with
  // respect to the point (xE,yE) given by the arithmetic mean.
  if ( verbose_flag ) { MSG("begin set_regn_geom_factors"<<endl<<flush) ; }
  assert(DIM==2) ;

  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {

    // get region vertices
    vector<int> vlist ;
    mesh.get_regn_vrtx( iR, vlist ) ;
    const int nRV = vlist.size() ;
    
    double xV[nRV], yV[nRV] ;
    double xE(0.), yE(0.) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = vlist[ilV] ;
      xV[ilV] = mesh.coords_V(iV,0) ;
      yV[ilV] = mesh.coords_V(iV,1) ;
      xE += xV[ilV] ;
      yE += yV[ilV] ;
    }
    xE /= double(nRV) ; 
    yE /= double(nRV) ; 
    
    // loop on the edges of the element
    double xR(0.), yR(0.), aR(0.) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      
      // quadrature stuff on the sub-triangle Ti (ilV,ilV+1,iR)
      // (vrtx quad rule exact for linear functions) 
      int il0 = ilV ;
      int il1 = (ilV+1)%nRV ;
      double xq[DIM+1] = { xE, xV[il0], xV[il1] } ;
      double yq[DIM+1] = { yE, yV[il0], yV[il1] } ;
      double aR_Ti = ( (xq[1]-xq[0])*(yq[2]-yq[0])-(xq[2]-xq[0])*(yq[1]-yq[0]) )/2. ;
      double wq    = aR_Ti/double(DIM+1) ; 
      
      // compute the contribution from Ti
      double xR_Ti(0.), yR_Ti(0.) ;
      for ( int iq=0; iq<DIM+1; ++iq ) {
	// accumulate contribution from (xq,yq)
	xR_Ti += wq * xq[iq] ;
	yR_Ti += wq * yq[iq] ;
      }
      xR += xR_Ti ;
      yR += yR_Ti ;
      aR += aR_Ti ;
    }

    // assign (xR,yR) and aR
    mesh.regn_coords[DIM*iR+0] = xR/aR ;
    mesh.regn_coords[DIM*iR+1] = yR/aR ;
    mesh.regn_area  [iR]       = aR ;
  }
  if ( verbose_flag ) { MSG("end-->set_regn_geom_factors"<<endl<<flush) ; }
}
void mesh2Dv_builder :: set_face_geom_factors() {
  if ( verbose_flag ) { MSG("begin set_face_geom_factors"<<endl<<flush) ; }
  assert(DIM==2) ;

  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {

    vector<int> vlist ;
    mesh.get_face_vrtx( iF, vlist ) ;
    assert( vlist.size()==DIM ) ;
    
    int iV0 = vlist[0] ;
    int iV1 = vlist[1] ;
    mesh.face_length[iF] = sqrt( pow( mesh.coords_V(iV1,0)-mesh.coords_V(iV0,0),2 ) +
				 pow( mesh.coords_V(iV1,1)-mesh.coords_V(iV0,1),2 ) ) ;

  } // end of --> for ( int iF=0; iF<nF; ++iF ) {...
  if ( verbose_flag ) { MSG("end-->set_face_geom_factors"<<endl<<flush) ; }
}

void mesh2Dv_builder :: change_bbox( double new_xmin, double new_ymin,      // min vertex
				     double new_xmax, double new_ymax  ) {  // max vertex
  assert( new_xmax>new_xmin ) ;
  assert( new_ymax>new_ymin ) ;
  assert( DIM==2 ) ;

  double new_max[DIM] = { new_xmax, new_ymax } ;
  double new_min[DIM] = { new_xmin, new_ymin } ;

  double xmin, ymin, xmax, ymax ;
  mesh.bbox( xmin, ymin, xmax, ymax ) ;

  double vol_scal = 1. ;
  for ( int s=0; s<DIM; ++s ) {
    vol_scal *= ( new_max[s]-new_min[s] ) / ( xmax-xmin ) ;
  }

  // set bounding box
  if ( !mesh.bb_status ) { mesh.eval_bbox() ; }

  double bb_den[DIM] ;
  for ( int s=0; s<DIM; ++s ) {
    bb_den[s] = mesh.bb_max[s]-mesh.bb_min[s] ;
  }

  // set new vertex coords
  for ( int iV=0; iV<mesh.nV; ++iV ) {
    for ( int s=0; s<DIM; ++s ) {
      mesh.vrtx_coords[DIM*iV+s] = ( ( mesh.vrtx_coords[DIM*iV+s]-mesh.bb_min[s] ) * new_max[s] + 
				     ( mesh.bb_max[s]-mesh.vrtx_coords[DIM*iV+s] ) * new_min[s] ) / bb_den[s] ;
    }
  } 
  
  // set new region coords
  for ( int iR=0; iR<mesh.nR; ++iR ) {
    for ( int s=0; s<DIM; ++s ) {
      mesh.regn_coords[DIM*iR+s] = ( ( mesh.regn_coords[DIM*iR+s]-mesh.bb_min[s] ) * new_max[s] + 
				     ( mesh.bb_max[s]-mesh.regn_coords[DIM*iR+s] ) * new_min[s] ) / bb_den[s] ;
    }
    mesh.regn_area[iR] *= vol_scal ;
  }
  
  // set face's geometric factors
  set_face_geom_factors() ;

  // reset bbox
  mesh.eval_bbox() ;
}
