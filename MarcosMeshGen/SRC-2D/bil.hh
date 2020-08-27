/*

  ONLY FOR QUADS: remap the unit square [0,1]x[0,1] into a quad with vertex coordinates
  (xvrt[0:3], yvrt[0:3])

 */
void build_bil_quad_mesh( mesh_2Dv & mesh, vector<double> & xvrt, vector<double> & yvrt, 
			 int nx=1, int ny=1, int mesh_flag=0, double wx=0., double wy=0. ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int cmap_flag = mesh_flag % 10 ;
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  int nV = (nx+1)*(ny+1) ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  double x, y ;
  double dx = 1./double(nx) ;
  double dy = 1./double(ny) ;
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      int iV = cmap.vrtx_index(i,j) ;
      double xh = double(i)*dx ;
      double yh = double(j)*dy ;
      cmap.bil_map_xy( xh, yh, xvrt, yvrt, x, y ) ; 
      xV[iV] = x ;
      yV[iV] = y ;
      fV[iV] = UNSET ;
      //cout << iV << " i=" << i << " j=" << j << " xV=" << xV[iV] << " yV=" << yV[iV] << endl ;
    }
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {
      regn_vlist.push_back( 4 ) ;
      regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
      regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
      regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
      regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
      fR.push_back( UNSET ) ;
    }
  }
  int nR = nx*ny ;
  regn_vlist.push_back( nR ) ;
  assert( nR == regn_vlist.size()/5 ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  char str_flag[10] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [10] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [10] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
