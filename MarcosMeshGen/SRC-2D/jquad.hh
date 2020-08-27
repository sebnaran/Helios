//-----------------------------------------------------------------------------------------
//      ---------------------------       ----------------------------
//      |     /|     |     /|     |       |     /|     |      |\     |
//      |    / |     |    / |     |       |    / |     |      | \    |
//      |   /  |     |   /  |     |       |   /  |     |      |  \   |
//      |  /   |     |  /   |     |       |  /   |     |      |   \  |
//      | /    |     | /    |     |       | /    |     |      |    \ |
//      |/     |     |/     |     |       |/     |     |      |     \|
//      ---------------------------       ----------------------------
//  
//      mesh_flag == (1)1x                   mesh_flag == (1)2x 
//
//-----------------------------------------------------------------------------------------
void build_jquad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
		       double wx=0., double wy=0. ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int cmap_flag = mesh_flag % 10 ;
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  int nV = 0 ;
  vector<int>    fV ;
  vector<double> xV, yV ;

  double dx = 1./double(nx) ;
  double dy = 1./double(ny) ;

  MatInt mask_vrtx( nx+1, ny+1 ) ;
  MatInt mask_edge( nx+1, ny+1 ) ;

  // counter 
  int iV = 0 ; 

  // nodes
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      //xV.push_back( dx * double(i) ) ; 
      //yV.push_back( dy * double(j) ) ; 
      double xij = dx * double(i) ;
      double yij = dy * double(j) ;
      xV.push_back( cmap.sxy( xij, yij ) ) ; 
      yV.push_back( cmap.txy( xij, yij ) ) ; 
      fV.push_back( UNSET ) ;
      mask_vrtx(i,j) = iV ;
      ++iV ;
      //cout << iV << " i=" << i << " j=" << j << " xV=" << xV[iV] << " yV=" << yV[iV] << endl ;
    }
  }
  // horizontal edges
  for ( int i=0; i<nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      //xV.push_back( dx * double(i) + dx/2. ) ; 
      //yV.push_back( dy * double(j)         ) ; 
      double xij = dx * double(i) + dx/2.  ;
      double yij = dy * double(j)          ;
      xV.push_back( cmap.sxy( xij, yij ) ) ; 
      yV.push_back( cmap.txy( xij, yij ) ) ; 
      fV.push_back( UNSET ) ;
      mask_edge(i,j) = iV ;
      ++iV ;
      //cout << iV << " i=" << i << " j=" << j << " xV=" << xV[iV] << " yV=" << yV[iV] << endl ;
    }
  }

  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {
      if ( mesh_flag/10==1 || 
	   ( mesh_flag/10==2 && i%2==0) ||
	   ( mesh_flag/10==3 && i%2==0 && j%2==0 ) ) { // 11?, by translation
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_vrtx(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_vrtx(i,j  ) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	regn_vlist.push_back( mask_vrtx(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( mask_edge(i,  j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j+1) ) ;
	regn_vlist.push_back( mask_edge(i,  j+1) ) ;
	fR.push_back( UNSET ) ;
      } else if ( ( mesh_flag/10==2 && i%2==1 ) || // 12? by reflection
		  ( mesh_flag/10==3 && i%2==1 && j%2==0 ) ) { 
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( mask_vrtx(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	regn_vlist.push_back( mask_vrtx(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_edge(i,  j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j  ) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_edge(i,  j+1) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j+1) ) ;
	fR.push_back( UNSET ) ;
      } else if ( mesh_flag/10==3 && i%2==0 && j%2==1 ) { // 15?, by double reflection
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_vrtx(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j) ) ;
	regn_vlist.push_back( mask_vrtx(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_edge(i,j  ) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	regn_vlist.push_back( mask_vrtx(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( mask_edge(i,  j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j+1) ) ;
	regn_vlist.push_back( mask_edge(i,  j+1) ) ;
	fR.push_back( UNSET ) ;
      } else if ( mesh_flag/10==3 && i%2==1 && j%2==1 ) { // 15?, by double reflection
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( mask_vrtx(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j) ) ;
	regn_vlist.push_back( mask_edge(i,j+1) ) ;
	regn_vlist.push_back( mask_vrtx(i,j+1) ) ;
	fR.push_back( UNSET ) ;
	
	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_edge(i,  j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j+1) ) ;
	fR.push_back( UNSET ) ;

	regn_vlist.push_back( 3 ) ;
	regn_vlist.push_back( mask_edge(i,  j  ) ) ;
	regn_vlist.push_back( mask_vrtx(i+1,j+1) ) ;
	regn_vlist.push_back( mask_edge(i  ,j+1) ) ;
	fR.push_back( UNSET ) ;
    } 
    }
  }
  int nR = 3*nx*ny ;
  regn_vlist.push_back( nR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  char str_flag[10] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [10] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [10] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("JQuad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
