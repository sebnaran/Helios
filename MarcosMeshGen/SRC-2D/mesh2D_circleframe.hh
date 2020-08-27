// mesh with circular strip (hardcoded)
void build_cquad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0, double wx=0., double wy=0. ) {
  MSGF("begin  build_cquad_mesh") ;
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int cmap_flag = mesh_flag % 10 ;
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  double xc = -0.5 ;
  double yc = +0.5 ;
  double dx = 1./double(nx) ;
  double dy = 1./double(ny) ;
  double xm = double(nx)*dx/2. ;
  // ---
  double x0 = 0.45 ;
  double x1 = 0.55 ;
  int    nV = (nx+1)*(ny+1) ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  for ( int i=0; i<=nx; ++i ) {
    // ---
    double si = dx*float(i) ; 
    double ri = sqrt( pow(si-xc,2) + pow(yc,2) ) ;
    // ---
    for ( int j=0; j<=ny; ++j ) {
      // ----
      double tj = dy*float(j) ; 
      // ----
      int    iV = cmap.vrtx_index(i,j) ;
      double th = atan( (tj-yc)/(si-xc) ) ;
      // ----
      double xi  = double(i)*dx ;
      double xcc = xc + ri*cos(th) ;
      double dxx = xcc - xi ;
      // ----
      double sc = 1. ;
      if      ( xi<x0 ) { sc = xi/x0           ; }
      else if ( xi>x1 ) { sc = (xi-1.)/(x1-1.) ; }
      if ( i==0 || i==nx || j==0 || j==ny ) { sc = 0. ; }
      // ----
      xV[iV] = xi + dxx * sc ;
      yV[iV] = yc + ri*sin(th) ;
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
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("CircleFrame-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
  // ----
  MSGF("end of build_cquad_mesh") ;
}
