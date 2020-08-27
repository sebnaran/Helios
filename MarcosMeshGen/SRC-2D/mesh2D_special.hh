#ifndef _MESH_2D_SPECIAL_HH
#define _MESH_2D_SPECIAL_HH

// 
// file for building some special meshes
// 
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// COORDINATE MAPS for REFLECTIONS of [0,1]x[0,1] domain
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*

  .    2
  ------------
  |          |
  |          |  
3 |          | 1
  |          |
  |          |
  ------------
  .    0
  
*/
class SideReflector {
private:
  const int    flag ;
  const double eps  ;
  bool is_equal( double a, double b ) const { return abs(a-b)<eps ; }

public:
  SideReflector( int _flag ) : flag(_flag), eps(1.e-12) {
    assert( 0<=flag && flag<=3 ) ;
  }
  ~SideReflector(){}

  bool is_on_reflection_boundary( double x, double y ) const {
    bool retval(false) ;
    switch( flag ) {
    case 0: retval = is_equal(y,0.) ; break ;
    case 1: retval = is_equal(x,1.) ; break ;
    case 2: retval = is_equal(y,1.) ; break ;
    case 3: retval = is_equal(x,0.) ; break ;
    }
    return retval ;
  }
  double new_X( double x ) const {
    double retval(0.) ;
    switch( flag ) {
    case 0: retval = x    ; break ;
    case 1: retval = 2.-x ; break ;
    case 2: retval = x    ; break ;
    case 3: retval = -x   ; break ;
    }
    return retval ;
  }
  double new_Y( double y ) const {
    double retval(0.) ;
    switch( flag ) {
    case 0: retval = -y   ; break ;
    case 1: retval = y    ; break ;
    case 2: retval = 2-y  ; break ;
    case 3: retval = y    ; break ;
    }
    return retval ;
  }
} ;
/*

 3            2 
  ------------
  |          |
  |          |  
  |          | 
  |          |
  |          |
  ------------
 0            1
  
 (i) reflection boundaries are boundaries incident to the reflection point
 Example:
 reflection on node 0, boundaries are 0-1 and 0-3 
 
 (ii) nodes on reflection boundaries are not included
*/
class PointReflector {
private:
  const int    flag ;
  const double eps  ;
  bool is_equal( double a, double b ) const { return abs(a-b)<eps ; }

public:
  PointReflector( int _flag ) : flag(_flag), eps(1.e-12) {
    assert( 0<=flag && flag<=3 ) ;
  }
  ~PointReflector(){}

  bool is_on_reflection_boundary( double x, double y ) const {
    bool retval(false) ;
    switch( flag ) {
    case 0: retval = is_equal(x,0.) || is_equal(y,0) ; break ;
    case 1: retval = is_equal(x,1.) || is_equal(y,0) ; break ;
    case 2: retval = is_equal(x,1.) || is_equal(y,1) ; break ;
    case 3: retval = is_equal(x,0.) || is_equal(y,1) ; break ;
    }
    return retval ;
  }
  double new_X( double x ) const {
    double retval(0.) ;
    switch( flag ) {
    case 0: retval = -x   ; break ;
    case 1: retval = 2.-x ; break ;
    case 2: retval = 2.-x ; break ;
    case 3: retval = -x   ; break ;
    }
    return retval ;
  }
  double new_Y( double y ) const {
    double retval(0.) ;
    switch( flag ) {
    case 0: retval = -y   ; break ;
    case 1: retval = -y   ; break ;
    case 2: retval = 2.-y ; break ;
    case 3: retval = 2.-y ; break ;
    }
    return retval ;
  }
} ;

class VectorTranslator {
private:
  const double vec_x, vec_y, eps ;
  bool is_equal( double a, double b ) const { return abs(a-b)<eps ; }
  
public:
  VectorTranslator( double _vec_x, double _vec_y ) : 
    vec_x(_vec_x), vec_y(_vec_y), eps(1.e-12) {}
  ~VectorTranslator(){}

  int is_on_side_boundary( double x, double y ) const {
    int retval = -1 ;
    if      ( is_equal(y,0.) ) { retval =  0 ; }
    else if ( is_equal(x,1.) ) { retval =  1 ; }
    else if ( is_equal(y,1.) ) { retval =  2 ; }
    else if ( is_equal(x,0.) ) { retval =  3 ; }
    else                       { retval = -1 ; }
    return retval ;
  }
  double new_X( double x ) const { return x + vec_x ; }
  double new_Y( double y ) const { return y + vec_y ; }
} ;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
//  LSHAPED-DOMAIN
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

namespace LSHAPED_DOMAIN {
  const double eps = 1.e-12 ;
  bool is_equal( double a, double b ) { return abs(a-b)<eps ; }
  bool is_point_on_boundary( double x, double y ) {
    bool is_horizontal_side = is_equal(y,-1.) || is_equal(y,1.) ||
      ( is_equal(y,0.) && 0.<= x && x<=1. ) ;
    bool is_vertical_side   = is_equal(x,-1.) || is_equal(x,1.) ||
      ( is_equal(x,0.) && -1.<= y && y<=1. ) ;
    return is_horizontal_side || is_vertical_side ;
  }
}

void build_input_for_Lshaped_domain_mesh( mesh_2Dv & mesh, 
					  vector<double> & xV, vector<double> & yV, vector<int> & fV,
					  vector<int>    & regn_vlist, vector<int> & fR ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  SideReflector  left_refl  (3) ;
  PointReflector origin_refl(0) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF A MESH COVERING an L-SHAPED DOMAIN
  //-----------------------------------------------------------------------------------------
  // first set of nodes
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    xV.push_back( mesh.coords_V(iV,0) ) ;
    yV.push_back( mesh.coords_V(iV,1) ) ;
    fV.push_back( mesh.get_fV(iV) ) ;
  }
  // nodes that are reflected by the left side
  int new_iV = mesh.n_vertex() ;
  int mask_left_nodes[mesh.n_vertex()] ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = left_refl.new_X( mesh.coords_V(iV,0) ) ;
    double yv = left_refl.new_Y( mesh.coords_V(iV,1) ) ;
    if ( !left_refl.is_on_reflection_boundary(xv,yv) ) {
      xV.push_back( xv ) ;
      yV.push_back( yv ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
      mask_left_nodes[iV] = new_iV ;
      ++new_iV ;
    } else {
      mask_left_nodes[iV] = iV ;
    }
  }
  // nodes that are reflected by the origin
  int mask_origin_nodes[mesh.n_vertex()] ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = origin_refl.new_X( mesh.coords_V(iV,0) ) ;
    double yv = origin_refl.new_Y( mesh.coords_V(iV,1) ) ;
    if ( !origin_refl.is_on_reflection_boundary(xv,yv) ) {
      xV.push_back( xv ) ;
      yV.push_back( yv ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
      mask_origin_nodes[iV] = new_iV ;
      ++new_iV ;
    } else {
      if ( abs(yv)<1.e-12 ) { // TOP side boundary
	mask_origin_nodes[iV] = mask_left_nodes[iV] ;
      } else if ( abs(xv)<1.e-12 ) { // RIGHT side boundary 
	xV.push_back( xv ) ;
	yV.push_back( yv ) ;
	fV.push_back( mesh.get_fV(iV) ) ;
	mask_origin_nodes[iV] = new_iV ;
	++new_iV ;
      } 
      else {
	assert(false) ;
      }
    }
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> loc_regn_vlist ;
    mesh.get_regn_vrtx( iR, loc_regn_vlist ) ;
    int nRV = loc_regn_vlist.size() ;
    // region in the first domain
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      regn_vlist.push_back( loc_regn_vlist[ilV] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
    // region in the left-reflected domain
    // set in counterclock-wise way
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=nRV-1; ilV>=0; --ilV ) {
      regn_vlist.push_back( mask_left_nodes[ loc_regn_vlist[ilV] ] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
    // region in the origin-reflected domain
    // set in counterclock-wise way
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      regn_vlist.push_back( mask_origin_nodes[ loc_regn_vlist[ilV] ] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
  }
  int nR = 3*mesh.n_region() ;
  regn_vlist.push_back( nR ) ;
}

void build_Lshaped_domain_mesh( mesh_2Dv & mesh ) {
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX/REGION DATA OF A MESH COVERING an L-SHAPED DOMAIN
  //-----------------------------------------------------------------------------------------
  vector<int>    fV ;
  vector<double> xV, yV ;
  vector<int> regn_vlist, fR ;
  build_input_for_Lshaped_domain_mesh( mesh, xV, yV, fV, regn_vlist, fR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  // build new mesh name from old one (standard flag 1xx)
  string new_mesh_name = string("L-Shaped-domain-")+mesh.get_mesh_name() ;
  // build new mesh
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set new mesh name
  mesh . set_mesh_name( new_mesh_name ) ;
}

void build_Lshaped_domain_mesh( int ilev, mesh_2Dv & mesh, Read_Mesh2D_Pars & rpars ) {
  //-----------------------------------------------------------------------------------------
  // build preliminar (regular) mesh
  //-----------------------------------------------------------------------------------------
  int mesh_flag = rpars.get_mesh_flag() ;
  int inpt_flag = abs(mesh_flag) ;
  int mf = inpt_flag-100*(inpt_flag/100) ;
  int nl = int( pow(2.,ilev) ) ; 
  int nx = rpars.get_mesh_nx() * nl ; 
  int ny = rpars.get_mesh_ny() * nl ; 
  switch( inpt_flag/100 ) {
  case 0: read_input_mesh( mesh, rpars ) ;              break ; // 0xx --> read from disk    
  case 1: build_quad_mesh ( mesh, nx, ny, mf, 0., 0. ) ; break ; // 1xx --> quads
  case 2: build_tria_mesh ( mesh, nx, ny, mf, 0., 0. ) ; break ; // 2xx --> triangles
  default: assert(false) ;
  }
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX/REGION DATA OF A MESH COVERING an L-SHAPED DOMAIN
  //-----------------------------------------------------------------------------------------
  vector<double> xV, yV ;
  vector<int>    fV, regn_vlist, fR ;
  //-----------------------------------------------------------------------------------------
  // build data for L-SHAPED DOMAIN MESH
  //-----------------------------------------------------------------------------------------
  build_input_for_Lshaped_domain_mesh( mesh, xV, yV, fV, regn_vlist, fR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  // build new mesh name from old one (standard flag 1xx)
  string new_mesh_name = string("L-Shaped-domain-")+mesh.get_mesh_name() ;
  // build new mesh
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set new mesh name
  mesh . set_mesh_name( new_mesh_name ) ;
  //-----------------------------------------------------------------------------------------
  // reset vertex coordinates
  //-----------------------------------------------------------------------------------------
  int cmap_flag = inpt_flag % 10 ;
  double wx = rpars.get_mesh_wx() ; 
  double wy = rpars.get_mesh_wy() ; 
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double   xv = mesh.coords_V( iV, 0 ) ;
    double   yv = mesh.coords_V( iV, 1 ) ;
    bool ok_remap = !LSHAPED_DOMAIN::is_point_on_boundary( xv, yv ) ;
    double new_xv = ok_remap ? cmap.sxy( xv, yv ) : xv ;
    double new_yv = ok_remap ? cmap.txy( xv, yv ) : yv ;
    mesh.set_vrtx_coords( iV, new_xv, new_yv ) ;
  }
  //-----------------------------------------------------------------------------------------
  // generation of dual mesh, if required
  //-----------------------------------------------------------------------------------------
  if ( mesh_flag<0 ) {
    dual_mesh_generation( mesh, mesh, rpars ) ;
  }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
// DOMAIN built by using four quadrants
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

namespace FOURQUADS_DOMAIN {
  const double eps = 1.e-12 ;
  bool is_equal( double a, double b ) { return abs(a-b)<eps ; }
  bool is_point_on_boundary( double x, double y ) {
    bool is_horizontal_side = is_equal(y,-1.) || is_equal(y,1.) ;
    bool is_vertical_side   = is_equal(x,-1.) || is_equal(x,1.) ;
    bool is_internal_side   = is_equal(x, 0.) || is_equal(y,0.) ;  
    return is_horizontal_side || is_vertical_side || is_internal_side ;
  }
}

void build_input_for_FourQuads_domain_mesh( mesh_2Dv & mesh, 
					    vector<double> & xV, vector<double> & yV, vector<int> & fV,
					    vector<int>    & regn_vlist, vector<int> & fR ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  SideReflector  left_refl  (3) ;
  SideReflector  bottom_refl(0) ;
  PointReflector origin_refl(0) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF A MESH COVERING an FOUR QUADRANTS-SHAPED DOMAIN
  //-----------------------------------------------------------------------------------------
  // first set of nodes
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    xV.push_back( mesh.coords_V(iV,0) ) ;
    yV.push_back( mesh.coords_V(iV,1) ) ;
    fV.push_back( mesh.get_fV(iV) ) ;
  }
  // nodes that are reflected by the left side
  int new_iV = mesh.n_vertex() ;
  int mask_left_nodes[mesh.n_vertex()] ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = left_refl.new_X( mesh.coords_V(iV,0) ) ;
    double yv = left_refl.new_Y( mesh.coords_V(iV,1) ) ;
    if ( !left_refl.is_on_reflection_boundary(xv,yv) ) {
      xV.push_back( xv ) ;
      yV.push_back( yv ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
      mask_left_nodes[iV] = new_iV ;
      ++new_iV ;
    } else {
      mask_left_nodes[iV] = iV ;
    }
  }
  // nodes that are reflected by the bottom side
  int mask_bottom_nodes[mesh.n_vertex()] ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = bottom_refl.new_X( mesh.coords_V(iV,0) ) ;
    double yv = bottom_refl.new_Y( mesh.coords_V(iV,1) ) ;
    if ( !bottom_refl.is_on_reflection_boundary(xv,yv) ) {
      xV.push_back( xv ) ;
      yV.push_back( yv ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
      mask_bottom_nodes[iV] = new_iV ;
      ++new_iV ;
    } else {
      mask_bottom_nodes[iV] = iV ;
    }
  }
  // nodes that are reflected by the origin
  int mask_origin_nodes[mesh.n_vertex()] ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = origin_refl.new_X( mesh.coords_V(iV,0) ) ;
    double yv = origin_refl.new_Y( mesh.coords_V(iV,1) ) ;
    if ( !origin_refl.is_on_reflection_boundary(xv,yv) ) {
      xV.push_back( xv ) ;
      yV.push_back( yv ) ;
      fV.push_back( mesh.get_fV(iV) ) ;
      mask_origin_nodes[iV] = new_iV ;
      ++new_iV ;
    } else {
      if ( abs(yv)<1.e-12 ) { // bottom side boundary
	mask_origin_nodes[iV] = mask_left_nodes[iV] ;
      } else if ( abs(xv)<1.e-12 ) { // left side boundary 
	mask_origin_nodes[iV] = mask_bottom_nodes[iV] ;
      } else {
	assert(false) ;
      }
    }
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> loc_regn_vlist ;
    mesh.get_regn_vrtx( iR, loc_regn_vlist ) ;
    int nRV = loc_regn_vlist.size() ;
    // region in the first domain
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      regn_vlist.push_back( loc_regn_vlist[ilV] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
    // region in the left-reflected domain
    // set in counterclock-wise way
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=nRV-1; ilV>=0; --ilV ) {
      regn_vlist.push_back( mask_left_nodes[ loc_regn_vlist[ilV] ] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
    // region in the bottom-reflected domain
    // set in counterclock-wise way
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=nRV-1; ilV>=0; --ilV ) {
      regn_vlist.push_back( mask_bottom_nodes[ loc_regn_vlist[ilV] ] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
    // region in the origin-reflected domain
    // set in counterclock-wise way
    regn_vlist.push_back( nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      regn_vlist.push_back( mask_origin_nodes[ loc_regn_vlist[ilV] ] ) ;
    }
    fR.push_back( mesh.get_fR(iR) ) ;
  }
  int nR = 4*mesh.n_region() ;
  regn_vlist.push_back( nR ) ;
}

void build_FourQuads_domain_mesh( mesh_2Dv & mesh ) {
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX/REGION DATA OF A MESH COVERING an FourQuads DOMAIN
  //-----------------------------------------------------------------------------------------
  vector<int>    fV ;
  vector<double> xV, yV ;
  vector<int> regn_vlist, fR ;
  //-----------------------------------------------------------------------------------------
  // build data for FourQuads DOMAIN MESH
  //-----------------------------------------------------------------------------------------
  build_input_for_FourQuads_domain_mesh( mesh, xV, yV, fV, regn_vlist, fR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  // build new mesh name from old one (standard flag 1xx)
  string new_mesh_name = string("FourQuads-domain-")+mesh.get_mesh_name() ;
  // build new mesh
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set new mesh name
  mesh . set_mesh_name( new_mesh_name ) ;
}

void build_FourQuads_domain_mesh( int ilev, mesh_2Dv & mesh, Read_Mesh2D_Pars & rpars ) {
  //-----------------------------------------------------------------------------------------
  // build preliminar (regular) mesh
  //-----------------------------------------------------------------------------------------
  int mesh_flag = rpars.get_mesh_flag() ;
  int inpt_flag = abs(mesh_flag) ;
  int mf = inpt_flag-100*(inpt_flag/100) ;
  int nl = int( pow(2.,ilev) ) ; 
  int nx = rpars.get_mesh_nx() * nl ; 
  int ny = rpars.get_mesh_ny() * nl ; 
  switch( inpt_flag/100 ) {
  case 1: build_quad_mesh ( mesh, nx, ny, mf, 0., 0. ) ; break ; // 1xx --> quads
  case 2: build_tria_mesh ( mesh, nx, ny, mf, 0., 0. ) ; break ; // 2xx --> triangles
  default: assert(false) ;
  }
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX/REGION DATA OF A MESH COVERING a Four Quad-Squared DOMAIN
  //-----------------------------------------------------------------------------------------
  vector<double> xV, yV ;
  vector<int>    fV, regn_vlist, fR ;
  //-----------------------------------------------------------------------------------------
  // build data for FourQuads DOMAIN MESH
  //-----------------------------------------------------------------------------------------
  build_input_for_FourQuads_domain_mesh( mesh, xV, yV, fV, regn_vlist, fR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  // build new mesh name from old one (standard flag 1xx)
  string new_mesh_name = string("FourQuads-domain-")+mesh.get_mesh_name() ;
  // build new mesh
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set new mesh name
  mesh . set_mesh_name( new_mesh_name ) ;
  //-----------------------------------------------------------------------------------------
  // reset mesh coordinates
  //-----------------------------------------------------------------------------------------
  int cmap_flag = inpt_flag % 10 ;
  double wx = rpars.get_mesh_wx() ; 
  double wy = rpars.get_mesh_wy() ; 
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xv = mesh.coords_V(iV,0) ;
    double yv = mesh.coords_V(iV,1) ;
    bool ok_remap = !FOURQUADS_DOMAIN::is_point_on_boundary( xv, yv ) ;
    double new_xv = ok_remap ? cmap.sxy( xv, yv ) : xv ;
    double new_yv = ok_remap ? cmap.txy( xv, yv ) : yv ;
    mesh.set_vrtx_coords( iV, new_xv, new_yv ) ;
  }
  // -- compute/set last geometric quantities
  mesh_builder.recompute_geom_factors() ;
  //-----------------------------------------------------------------------------------------
  // generation of dual mesh, if required
  //-----------------------------------------------------------------------------------------
  if ( mesh_flag<0 ) {
    dual_mesh_generation( mesh, mesh, rpars ) ;
  }
  //-----------------------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
// SINGLE CELL MESH, REGULAR POLYGON WITH nvrtx VERTICES
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void build_regular_polygon_mesh( mesh_2Dv & mesh, int nvrtx=3 ) {
  assert( nvrtx>=3 ) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF REGULAR POLYGONAL MESH
  //-----------------------------------------------------------------------------------------
  double pig = acos(-1.) ;
  double cff = 2.*pig/double(nvrtx) ;
  // ---
  int nV = nvrtx ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  for ( int iV=0; iV<nV; ++iV ) {
    xV[iV] = cos( cff*double(iV) ) ;
    yV[iV] = sin( cff*double(iV) ) ;
    fV[iV] = UNSET ;
    printf("%i --> (%14.7e,%14.7e)\n",iV,xV[iV],yV[iV]) ;
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF REGULAR POLYGONAL MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  regn_vlist.push_back( nV ) ;  
  for ( int iV=0; iV<nV; ++iV ) {
    regn_vlist.push_back( iV ) ;
  }
  fR.push_back( UNSET ) ;
  int nR = 1 ;
  regn_vlist.push_back( nR ) ;
  assert( nR == regn_vlist.size()/(1+nV) ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR REGULAR POLYGONAL MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  string mesh_name = string("Regular Polygonal Mesh") ;
  mesh . set_mesh_name( mesh_name ) ;
}

#endif // end of _MESH_2D_SPECIAL_HH
