//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// QUADRILATERAL-BASED MESH (standard flag: 1xx)
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class VertexIndex {
private:
  int i0, j0, i1, j1 ;
  int nx, ny, nV ;
  vector<int> vrtx_index ;
  vector<int> vrtx_flag  ;
  int idx( int i, int j ) const {
    return i + (nx+1) * j ;
  }
public:
  VertexIndex( int _nx, int _ny ) : nx(_nx), ny(_ny) {
    i0 = nx/2-1 ;
    i1 = nx/2+1 ;
    j0 = ny/2-1 ;
    j1 = ny/2+1 ;
    // set vrtx_index stuff
    vrtx_index.resize( (nx+1)*(ny+1) ) ;
    int iV = 0 ;
    for ( int i=0; i<=nx; ++i ) {
      for ( int j=0; j<=ny; ++j ) {
	if ( !( i0<i && i<i1 && j0<j && j<j1 ) ) {
	  vrtx_index[ idx(i,j) ] = iV++ ;
	} else {
	  vrtx_index[ idx(i,j) ] = -1 ;
	}
      }
    }
    nV = iV ;
    // set vrtx_flag
    vrtx_flag .resize( (nx+1)*(ny+1) ) ;
    for ( int i=0; i<=nx; ++i ) {
      for ( int j=0; j<=ny; ++j ) {
	vrtx_flag[ idx(i,j) ] = 0 ;
	if ( i==0 || i==nx || j==0 || j==ny ) {
	  vrtx_flag[ idx(i,j) ] = 1 ;
	}
	if ( ( i==i0 || i==i1 ) && ( j0<=j && j<=j1) ) {
	  vrtx_flag[ idx(i,j) ] = 2 ;
	}
	if ( ( j==j0 || j==j1 ) && ( i0<=i && i<=i1) ) {
	  vrtx_flag[ idx(i,j) ] = 2 ;
	}
      }
    }
  }
  ~VertexIndex(){}
  int operator() ( int i, int j ) {
    assert( 0<=i && i<=nx ) ;
    assert( 0<=j && j<=ny ) ;
    return vrtx_index[ idx(i,j) ] ;
  }
  int  n_vertex() { return nV ; }
  bool is_vrtx_in_hole( int i, int j ) { return i0<i  && i<i1 && j0<j  && j<j1 ; }
  bool is_regn_in_hole( int i, int j ) { return i0<=i && i<i1 && j0<=j && j<j1 ; }
  int  get_flag( int i, int j ) const {
    return vrtx_flag[ idx(i,j) ] ;
  }
  void prt_debug() {
    PRT(nV) ;
    PRT((nx+1)*(ny+1)) ;
    PRT(vrtx_index.size()) ;
    LINE(--) ;
    VAL(i0) ; VAL(i1) ; VAL(j0) ; PRT(j1) ;
    LINE(--) ;
    for ( int i=0; i<=nx; ++i ) {
      for ( int j=0; j<=ny; ++j ) {
	cout << " (" << i << "," << j << "), iV = " << vrtx_index[ idx(i,j) ] << ", fV = " << vrtx_flag[ idx(i,j) ] << endl ;
      }
    }
    LINE(--) ;
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j ) {
	cout << " (" << i << "," << j << "), region = " << is_regn_in_hole(i,j) << endl ;
      }
    }
  }
} ;
void build_hquad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
		       double wx=0., double wy=0. ) {
  
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int cmap_flag = mesh_flag % 10 ;
  VertexIndex vrtx_index(nx,ny) ;
  //vrtx_index.prt_debug() ;
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  int nV = vrtx_index.n_vertex() ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      int iV = vrtx_index(i,j) ;
      assert( iV==-1 || ( 0<=iV && iV<nV ) ) ;
      if ( 0<=iV && iV<nV ) { 
	xV[iV] = cmap.sxy(i,j) ;
	yV[iV] = cmap.txy(i,j) ;
	fV[iV] = vrtx_index.get_flag( i, j ) ;
      }
      //cout << iV << " i=" << i << " j=" << j << " xV=" << xV[iV] << " yV=" << yV[iV] << endl ;
    }
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  int nR = 0 ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {
      if ( !vrtx_index.is_regn_in_hole(i,j) ) {
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( vrtx_index(i,    j) ) ;
	regn_vlist.push_back( vrtx_index(i+1,  j) ) ;
	regn_vlist.push_back( vrtx_index(i+1,j+1) ) ;
	regn_vlist.push_back( vrtx_index(i,  j+1) ) ;
	fR.push_back( UNSET ) ;
	++nR ;
      }
    }
  }
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
  string mesh_name = string("Hole-Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
