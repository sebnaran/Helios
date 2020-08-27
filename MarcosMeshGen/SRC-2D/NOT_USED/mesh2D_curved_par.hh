#ifndef _MESH2DC_HH
#define _MESH2DC_HH

class mesh_2Dc ;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// LOGICALLY RECTANGULAR (CURVED) MAPPINGS
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class CurvedMap {
private:
  static const int HORIZONTAL_LINE = 0 ;
  static const int VERTICAL_LINE   = 1 ;

private:
  // input members
  int nx    ; // number of partitions along X
  int ny    ; // number of partitions along Y
  int mflag ; // mesh flag to select the mapping function
  double wi ; // mapping parameter along X
  double wj ; // mapping parameter along Y

  // local vars
  bool   status ;
  double pi, dx, dy ;
  double cx, cy ;

  // aux methods
  inline double xy_rand( double a=-1., double b=1. ) const {
    return a + (b-a) * double(random())/double(LONG_MAX) ;
  }
  inline bool is_bnd_vertex( int i, int j ) const {
    return i==0 || i==nx || j==0 || j==ny ;
  }
  inline bool is_i_on_boundary( int i ) const {
    return i==0 || i==nx ;
  }
  inline bool is_j_on_boundary( int j ) const {
    return j==0 || j==ny ;
  }

  // -------------------
  // map-0, identity map
  // -------------------
  inline double s0( int i, int j ) const { return dx * double(i) ; }
  inline double t0( int i, int j ) const { return dy * double(j) ; }
  // ---
  inline double s0  ( double x, double y ) const { return x  ; }
  inline double s0_x( double x, double y ) const { return 1. ; }
  inline double s0_y( double x, double y ) const { return 0. ; }
  // ---
  inline double t0  ( double x, double y ) const { return y  ; }
  inline double t0_x( double x, double y ) const { return 0. ; }
  inline double t0_y( double x, double y ) const { return 1. ; }

  // --------------------------------------
  // map-1, sin(..)sin(..)-Los Alamos map
  // --------------------------------------
  inline double s1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return sij + wi * sin(2.*pi*sij)*sin(2.*pi*tij) ; 
  }
  inline double t1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return tij + wj * sin(2.*pi*sij)*sin(2.*pi*tij) ; 
  }
  // 
  inline double s1( double x, double y ) const { 
    return x  + wi *       sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double s1_x( double x, double y ) const { 
    return 1. + wi * 2.*pi*cos(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double s1_y( double x, double y ) const { 
    return 0. + wi * 2.*pi*sin(2.*pi*x)*cos(2.*pi*y) ; 
  }

  inline double t1( double x, double y ) const { 
    return y  + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1_x( double x, double y ) const { 
    return 0. + wj * 2.*pi*cos(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1_y( double x, double y ) const { 
    return 1. + wj * 2.*pi*sin(2.*pi*x)*cos(2.*pi*y) ; 
  }

  // --------------------------------------
  // map-2, sin(..)sin(..)-Florence Hubert map
  // --------------------------------------
  inline double s2( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return sij ;
  }
  inline double t2( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return tij + wj * sin(wi*pi*sij) * tij*(1.-tij)*4. ;
  }
  // 
  inline double s2( double x, double y ) const { 
    return x  ;
  }
  inline double s2_x( double x, double y ) const { 
    return 1. ;
  }
  inline double s2_y( double x, double y ) const { 
    return 0. ;
  }

  inline double t2( double x, double y ) const { 
    return y  + wj * sin(wi*pi*x) *  y*(1.-y)*4. ;
  }
  inline double t2_x( double x, double y ) const { 
    return 0. + wj * wi*pi * cos(wi*pi*x) *  y*(1.-y)*4. ;
  }
  inline double t2_y( double x, double y ) const { 
    return 1. + wj * sin(wi*pi*x) * (1.-2.*y)*4. ;
  }
  
public:
  // (default is 1 cell, with identity map)
  CurvedMap( int _nx=1, int _ny=1, int _mflag=0, double _wi=0., double _wj=0. ) : 
    status(false) {
    setup(_nx,_ny,_mflag,_wi,_wj) ;
  }
  ~CurvedMap() {}
  inline void setup( int _nx, int _ny, int _mflag, double _wi, double _wj ) {
    assert( _nx>0 && _ny>0 ) ;
    // --- set status: if status==false, CurvedMap object will be unable to respond to 
    //                 mesh_2Dc queries
    status = true ;
    // --- store input data
    nx = _nx ;
    ny = _ny ;
    mflag = _mflag ;
    wi = _wi ;
    wj = _wj ;
    //--- set internal attributes
    dx = 1./double(nx) ;
    dy = 1./double(ny) ;
    pi = acos(-1.) ;
    // --- aux
    cx = 2. ;
    cy = wj ;
  }

  // return vrtx index
  inline int vrtx_index( int i, int j ) const {
    return i + (nx+1) * j ;
  }
  // MAPPINGS:
  // (i,j) --> ( csi_ij, eta_ij )
  inline double sxy( int i, int j ) const {
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=s0(i,j) ; break ;
    case 1: retval=s1(i,j) ; break ;
    case 2: retval=s2(i,j) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy( int i, int j ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=t0(i,j) ; break ;
    case 1: retval=t1(i,j) ; break ;
    case 2: retval=t2(i,j) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  // (x,y) --> csi(x,y), and derivatives
  inline double sxy( double x, double y ) const {
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=s0(x,y) ; break ;
    case 1: retval=s1(x,y) ; break ;
    case 2: retval=s2(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double sxy_x( double x, double y ) const {
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=s0_x(x,y) ; break ;
    case 1: retval=s1_x(x,y) ; break ;
    case 2: retval=s2_x(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double sxy_y( double x, double y ) const {
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=s0_y(x,y) ; break ;
    case 1: retval=s1_y(x,y) ; break ;
    case 2: retval=s2_y(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  // (x,y) --> eta(x,y), and derivatives
  inline double txy( double x, double y ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=t0(x,y) ; break ;
    case 1: retval=t1(x,y) ; break ;
    case 2: retval=t2(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy_x( double x, double y ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=t0_x(x,y) ; break ;
    case 1: retval=t1_x(x,y) ; break ;
    case 2: retval=t2_x(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy_y( double x, double y ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=t0_y(x,y) ; break ;
    case 1: retval=t1_y(x,y) ; break ;
    case 2: retval=t2_y(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }

public:
  // ---------------------------------------------------------------------------------------
  // MAPPINGS: PUBLIC INTERFACE WITH mesh_2Dc OBJECTS
  // ---------------------------------------------------------------------------------------
  // GRIDLINES, they are built using the previous mappings
  // to remap a logically rectangular grid into physical domain
  // grid lines are identified by two data: 
  // gl = 0, 1     --> grid line flag, horizontal, vertical
  // gv = <double> --> grid line position
  inline double xs( double s, int gl, double gv ) const {
    assert( gl==HORIZONTAL_LINE || gl==VERTICAL_LINE ) ;
    assert( status ) ; 
    double retval = 0. ;
    switch( gl ) {
    case HORIZONTAL_LINE : retval = sxy( s, gv ) ; break ;
    case VERTICAL_LINE   : retval = sxy( gv, s ) ; break ;
    }
    return retval ;
  }
  inline double ys( double s, int gl, double gv ) const {
    assert( gl==HORIZONTAL_LINE || gl==VERTICAL_LINE ) ;
    assert( status ) ; 
    double retval = 0. ;
    switch( gl ) {
    case HORIZONTAL_LINE : retval = txy( s, gv ) ; break ;
    case VERTICAL_LINE   : retval = txy( gv, s ) ; break ;
    }
    return retval ;
  }
  // derivatives along the gridlines
#if 1
  inline double xs_s( double s, int gl, double gv ) const {
    assert( gl==HORIZONTAL_LINE || gl==VERTICAL_LINE ) ;
    assert( status ) ; 
    double retval = 0. ;
    switch( gl ) {
    case HORIZONTAL_LINE : retval = sxy_x( s, gv ) ; break ;
    case VERTICAL_LINE   : retval = sxy_y( gv, s ) ; break ;
    }
    return retval ;
  }
  inline double ys_s( double s, int gl, double gv ) const {
    assert( gl==HORIZONTAL_LINE || gl==VERTICAL_LINE ) ;
    assert( status ) ; 
    double retval = 0. ;
    switch( gl ) {
    case HORIZONTAL_LINE : retval = txy_x( s, gv ) ; break ;
    case VERTICAL_LINE   : retval = txy_y( gv, s ) ; break ;
    }
    return retval ;
  }
#else // introduce an additional error
  inline double xs_s( double s, int gl, double gv ) const {
    double eps = 1.e-8 ;
    return ( xs(s+eps,gl,gv) - xs(s-eps,gl,gv) )/( 2.*eps ) ;
  }
  inline double ys_s( double s, int gl, double gv ) const {
    double eps = 1.e-8 ;
    return ( ys(s+eps,gl,gv) - ys(s-eps,gl,gv) )/( 2.*eps ) ;
  }
#endif
  // ---------------------------------------------------------------------------------------
  // end of public interface with mesh_2Dc objects
  // ---------------------------------------------------------------------------------------

  inline double get_dx() { return dx ; }
  inline double get_dy() { return dy ; }
} ;

#include "coax.hh"

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
// A face of a logically rectangular  cell.
// It needs the gridline id and (starting,ending) points (s0,s1)
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class curved_face {
private:
  int    grid_line ; // type of grid lines
  double grid_val  ; // grid line position
  double s0, s1 ;    // global parameterization (from s0 to s1)
  bool   status ;    // status
public:
  curved_face() : grid_line(-1), s0(0.), s1(0.), status(false) {}
  curved_face( int _grid_line, double _s0, double _s1 ) : 
    grid_line(_grid_line), s0(_s0), s1(_s1), status(true) {}
  ~curved_face() {}
  
  bool is_curved_face() { return status ; } 
  
  void setup( int _grid_line, double _grid_val, double _s0, double _s1 ) {
    grid_line = _grid_line ;
    grid_val  = _grid_val  ;
    s0        = _s0 ;
    s1        = _s1 ;
    status    = true ;
  }
  inline double get_s0()        { return s0 ; }
  inline double get_s1()        { return s1 ; }
  inline double get_grid_val () { return grid_val  ; }
  inline int    get_grid_line() { return grid_line ; }
} ;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//
// MESH MANAGER
//
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class mesh_2Dc : public mesh_2Dv {

private:
  // aux geometrical quantities
  typedef TDataMat<double> dataMat ;
  // ---
  dataMat   cvd_R_coords ;
  dataMat   cvd_R_area   ;
  // ---
  dataMat   cvd_F_length ;
  dataMat   cvd_F_nor_x ;
  dataMat   cvd_F_nor_y ;

  // possible curved parameterization
  int       cvd_param_flag ;
  CurvedMap cmap_0 ;
  CoaxMap   cmap_1 ;

public:
  void setup_CurvedMap( int _nx, int _ny, int _mflag, double _wi, double _wj ) {
    MSGF("working with CurvedMap") ;
    cvd_param_flag = 0 ;
    cmap_0.setup(_nx,_ny,_mflag,_wi,_wj) ;
  }
  void setup_CoaxMap( CoaxStructure & coax ) {
    MSGF("working with CoaxMap") ;
    cvd_param_flag = 1 ;
    cmap_1.setup( coax ) ;
  }

private:
  double xs  ( double s, int gl, double gv ) {
    double retval ;
    switch( cvd_param_flag ) {
    case 0 : retval = cmap_0.xs(s,gl,gv) ; break ;
    case 1 : retval = cmap_1.xs(s,gl,gv) ; break ;
    default : assert(false) ;
    }
    return retval ;
  }
  double ys  ( double s, int gl, double gv ) {
    double retval ;
    switch( cvd_param_flag ) {
    case 0 : retval = cmap_0.ys(s,gl,gv) ; break ;
    case 1 : retval = cmap_1.ys(s,gl,gv) ; break ;
    default : assert(false) ;
    }
    return retval ;
  }
  double xs_s( double s, int gl, double gv ) {
    double retval ;
    switch( cvd_param_flag ) {
    case 0 : retval = cmap_0.xs_s(s,gl,gv) ; break ;
    case 1 : retval = cmap_1.xs_s(s,gl,gv) ; break ;
    default : assert(false) ;
    }
    return retval ;
  }
  double ys_s( double s, int gl, double gv ) {
    double retval ;
    switch( cvd_param_flag ) {
    case 0 : retval = cmap_0.ys_s(s,gl,gv) ; break ;
    case 1 : retval = cmap_1.ys_s(s,gl,gv) ; break ;
    default : assert(false) ;
    }
    return retval ;
  }

private:
  bool cvd_status ;
  vector<curved_face> cvd_flist ;

  void setup_face_geometry() ;
  void setup_region_geometry() ;
  void setup_region_geometry( int iR ) ;

  // setup geometry of curved faces
  void check_region_geometry() ;
  void check_region_geometry( int iR ) ;
  
private:
  // thes quantities are for a global check of consistency
  double max_diff_area, max_diff_bary_X, max_diff_bary_Y ;

public:
  mesh_2Dc() : cvd_status(false), cvd_param_flag(-1) {}
  ~mesh_2Dc() {}

  // used in main to init curved faces
  void init_curve() ;
  void set_cvd_status( bool new_cvd_status ) ;

  // SHARED with mesh_2Dv
  // return true of iF is a curved face
  // total number of curved faces of the mesh and of region iR
  inline virtual int  n_curved_faces() ;
  inline virtual int  n_curved_faces( int iR ) ;
  inline virtual bool is_curved_face( int iF ) ;
  inline virtual bool is_curved_regn( int iR ) ;

  // get the list of curved faces of region iR
  void get_regn_cvd_face( int iR, vector<int> & cvd_flist ) ;

  // local shape of a curved face
  double face_pnts( double t, int iF, int idim ) ;
  double face_derv( double t, int iF, int idim ) ;

  // re-parameterization
  double face_alpha( double t, int iF ) ;
  double face_delta( double t, int iF ) ;

  // barycenter of region R and face iF (possibly curved)
  double cvd_coords_R( int iR, int k ) ;
  double cvd_coords_F( int iF, int k ) ;
  
  // area of region R (possibly with curved faces)
  double get_cvd_regn_measure( int iR ) ;

  // length of face F (possibly curved)
  double get_cvd_face_measure( int iF ) ;
  double get_cvd_nor( int iF, int s ) ;
  
  // initialization methods
  void setup_cvd_flist() {
    cvd_flist.resize( n_face() ) ;
  }
  void setup_cvd_flist( int iF, int grid_line, double grid_val, double s0, double s1 ) {
    assert( 0<=iF & iF<n_face() ) ;
    cvd_flist[iF].setup( grid_line, grid_val, s0, s1 ) ;
  }

  // DEBUG
  void check_mesh() {
    int nF = n_face() ;
    for ( int iF=0; iF<nF; ++iF ) {
      LINE(--) ;
      VAL(iF) ; PRT( is_curved_face(iF) ) ;
      if ( is_curved_face(iF) ) {
	double xA = face_pnts( 0., iF, 0 ) ;
	double yA = face_pnts( 0., iF, 1 ) ;
	double xB = face_pnts( 1., iF, 0 ) ;
	double yB = face_pnts( 1., iF, 1 ) ;
	VAL(xA) ; PRT(yA) ; 
	VAL(xB) ; PRT(yB) ; 
	double xm = face_pnts( 0.5, iF, 0 ) ;
	double ym = face_pnts( 0.5, iF, 1 ) ;
	VAL(xm) ; PRT(ym) ; 
	double xs0 = face_derv( 0., iF, 0 ) ;
	double ys0 = face_derv( 0., iF, 1 ) ;
	double xs1 = face_derv( 1., iF, 0 ) ;
	double ys1 = face_derv( 1., iF, 1 ) ;
	VAL(xs0) ; PRT(ys0) ; 
	VAL(xs1) ; PRT(ys1) ; 
      }
    }
  }

  void print_region_data() ;
  void print_region_data( int iR ) ;
} ;

// used in main to init curved faces
void mesh_2Dc :: init_curve() {
  if ( !cvd_status ) { 
    cvd_status = true ;
    setup_cvd_flist() ;
  }
  setup_region_geometry() ;
  check_region_geometry() ;
}

void mesh_2Dc ::  set_cvd_status( bool new_cvd_status ) {
  cvd_status = new_cvd_status ;
}

// every face has two extrema, so TWO points is the minimum
bool mesh_2Dc :: is_curved_face( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  return cvd_flist[iF].is_curved_face() ;
}

// check whther region iR has at least one curved face
bool mesh_2Dc :: is_curved_regn( int iR ) {
  assert( 0<=iR && iR<nR ) ;
  vector<int> regn_flist ;
  get_regn_face( iR, regn_flist ) ;
  int nRF = regn_flist.size() ;
  bool retval = false ;
  for ( int ilF=0; ilF<nRF && !retval; ++ilF ) {
    int iF = regn_flist[ilF] ;
    retval |= cvd_flist[iF].is_curved_face() ;
  }
  return retval ;
}

// total number of curved face of the mesh
int mesh_2Dc :: n_curved_faces() {
  int ncf = 0 ; // #curved faces
  for ( int iF=0; iF<n_face(); ++iF ) {
    if ( is_curved_face(iF) ) { ncf++ ; }
  }
  return ncf ;
}

// total number of curved faces of region iR
int mesh_2Dc :: n_curved_faces( int iR ) {
  assert( 0<=iR && iR<n_region() ) ;
  int ncf = 0 ; // #curved faces
  vector<int> regn_flist ;
  get_regn_face( iR, regn_flist ) ;
  for ( int ilF=0; ilF<regn_flist.size(); ++ilF ) {
    if ( is_curved_face( regn_flist[ilF] ) ) { ncf++ ; }
  }
  return ncf ;
}

// get the list of curved faces of region iR
void mesh_2Dc :: get_regn_cvd_face( int iR, vector<int> & cvd_flist ) {
  assert( 0<=iR && iR<n_region() ) ;
  cvd_flist.resize(0) ;
  vector<int> regn_flist ;
  get_regn_face( iR, regn_flist ) ;
  for ( int ilF=0; ilF<regn_flist.size(); ++ilF ) {
    int iF = regn_flist[ilF] ;
    if ( is_curved_face( iF ) ) { cvd_flist.push_back( iF ) ; }
  }
}

// parameterization of the curved face F through X(s)
// s in [s0,s1] --> X(s)
double mesh_2Dc :: face_pnts( double z, int iF, int idim ) {
  //assert( 0<=t    && t<=1        ) ;
  assert( 0<=iF   && iF<n_face() ) ;
  assert( 0<=idim && idim<DIM    ) ;
  double retval = 0. ;
  if ( is_curved_face(iF) ) {
    int    gl = cvd_flist[iF].get_grid_line() ;
    double gv = cvd_flist[iF].get_grid_val () ;
    double s0 = cvd_flist[iF].get_s0() ;
    double s1 = cvd_flist[iF].get_s1() ;
    double s  = s0 * (1.-z) + s1 * z ;
    if ( idim==0 ) { retval = xs( s, gl, gv ) ; }
    else           { retval = ys( s, gl, gv ) ; }
  } else { // to treat the case of a planar face
    int iV0 = face_vrtx( iF, 0 ) ;
    int iV1 = face_vrtx( iF, 1 ) ;
    retval = coords_V( iV0, idim ) * (1.-z) + coords_V( iV1, idim ) * z ;
  }
  return retval ;
}

// X'(s) \in F, s \in [s0,s1]
double mesh_2Dc :: face_derv( double z, int iF, int idim ) {
  //assert( 0<=t    && t<=1        ) ;
  assert( 0<=iF   && iF<n_face() ) ;
  assert( 0<=idim && idim<DIM    ) ;
  double retval = 0. ;
  if ( is_curved_face(iF) ) {
    int    gl    = cvd_flist[iF].get_grid_line() ;
    double gv    = cvd_flist[iF].get_grid_val () ;
    double s0    = cvd_flist[iF].get_s0() ;
    double s1    = cvd_flist[iF].get_s1() ;
    double s     = s0 * (1.-z) + s1 * z ;
    double ds_dz = s1-s0 ;
    if ( idim==0 ) { retval = ds_dz * xs_s( s, gl, gv ) ; }
    else           { retval = ds_dz * ys_s( s, gl, gv ) ; }
  } else { // to treat the case of a planar face
    int iV0 = face_vrtx( iF, 0 ) ;
    int iV1 = face_vrtx( iF, 1 ) ;
    retval = coords_V( iV1, idim ) - coords_V( iV0, idim ) ;
  }
  return retval ;
}

// distance of the projection of X(t) onto the face from the first point
double mesh_2Dc :: face_alpha( double t, int iF ) {
  assert( 0<=t  && t<=1        ) ;
  assert( 0<=iF && iF<n_face() ) ;
  double tx = get_tng( iF, 0 ) ;
  double ty = get_tng( iF, 1 ) ;
  double lF = get_face_measure( iF ) ;
  double px = face_pnts(t,iF,0)-face_pnts(0.,iF,0) ;
  double py = face_pnts(t,iF,1)-face_pnts(0.,iF,1) ;
  return ( tx*px + ty*py )/lF ;
}

// distance of X(t) from the planar face
double mesh_2Dc :: face_delta( double t, int iF ) {
  assert( 0<=t  && t<=1        ) ;
  assert( 0<=iF && iF<n_face() ) ;
  assert( is_curved_face(iF)     ) ;
  double nx = get_nor( iF, 0 ) ;
  double ny = get_nor( iF, 1 ) ;
  double lF = get_face_measure( iF ) ;
  double px = face_pnts(t,iF,0)-face_pnts(0.,iF,0) ;
  double py = face_pnts(t,iF,1)-face_pnts(0.,iF,1) ;
  return ( nx*px + ny*py )/lF ;
}

// barycenter of region R (possibly with curved faces)
double mesh_2Dc :: cvd_coords_R( int iR, int k ) {
  assert( 0<=iR && iR<nR ) ;
  assert( 0<=k  && k<DIM ) ;
  return cvd_R_coords( iR, k ) ;
}

// midpoint of (possibly curved) face iF
double mesh_2Dc :: cvd_coords_F( int iF, int k ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=k  && k<DIM ) ;
  double retval = 0. ;
  return face_pnts( 0.5, iF, k ) ;
}

// area of region R  (possibly with curved faces)
double mesh_2Dc :: get_cvd_regn_measure( int iR ) {
  assert( 0<=iR && iR<nR ) ;
  assert( R_area.size()==nR && R_area.size_loc(iR)==1 ) ;
  return cvd_R_area(iR,0) ;
}

// length of face F  (possibly curved)
double mesh_2Dc :: get_cvd_face_measure( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  return cvd_F_length(iF,0) ;
}
double mesh_2Dc :: get_cvd_nor( int iF, int s ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=s  && s<DIM ) ;
  return s==0 ? cvd_F_nor_x(iF,0) : cvd_F_nor_y(iF,0) ;
}

void mesh_2Dc :: setup_region_geometry() {
  MSG("begin setup_curved_regn_geom"<<endl<<flush) ;
  cvd_F_length.setup(nF,1) ;
  cvd_F_nor_x.setup (nF,1) ;
  cvd_F_nor_y.setup (nF,1) ;
  setup_face_geometry() ;
  // ---
  cvd_R_area  .setup(nR,1) ;
  cvd_R_coords.setup(nR,DIM) ;
  for ( int iR=0; iR<nR; ++iR ) {
    setup_region_geometry( iR ) ;
  }
  MSG("end-->setup_curved_regn_geom"<<endl<<flush) ;
}

void mesh_2Dc :: setup_face_geometry() {
  //Gauss4 quad ;
  GaussLegendre8 quad ;
  const int nq = quad.get_nq() ;
  double sq[nq], wq[nq] ;
  quad.get_quadrule( nq, sq, wq ) ;
  // loop on mesh faces
  for ( int iF=0; iF<nF; ++iF ) {

    if ( is_curved_face(iF) ) {

      // HIGHER_ORDER approximation
      // lenF <-- +int_{[0,1]} sqrt( x'(t)^2+y'(t)^2 ) dt
      double nxF(0.), nyF(0.), lF(0.) ;


      for ( int iq=0; iq<nq; ++iq ) {
	double len = sqrt( pow( face_derv( sq[iq], iF, 0 ), 2) +
			   pow( face_derv( sq[iq], iF, 1 ), 2 ) ) ;
	lF  += wq[iq] * len ;
	nxF += wq[iq] * face_derv( sq[iq], iF, 1 )  /len ;
	nyF -= wq[iq] * face_derv( sq[iq], iF, 0 )  /len ;
      }
      
      cvd_F_length(iF,0) = lF ;
      cvd_F_nor_x (iF,0) = nxF / lF ;
      cvd_F_nor_y (iF,0) = nyF / lF ;

    } else {

      // a flat face has the length already calculated by mesh2Dv
      cvd_F_length(iF,0) = get_face_measure(iF) ;
      cvd_F_nor_x (iF,0) = get_nor(iF,0) ;
      cvd_F_nor_y (iF,0) = get_nor(iF,1) ;
	
    }
  }
}

void mesh_2Dc :: setup_region_geometry( int iR ) {
  assert( 0<=iR && iR<n_region() ) ;

  vector<int> regn_flist ;
  get_regn_face( iR, regn_flist ) ;
  const int nRF = n_regn_face( iR ) ;
  const int nRV = nRF ;
  const int nFV = DIM ;

  //Gauss4 quad ;
  GaussLegendre8 quad ;
  const int nq = quad.get_nq() ;
  double sq[nq], wq[nq] ;
  quad.get_quadrule( nq, sq, wq ) ;

  double area_X(0.), area_Y(0.) ;
  double cntr_X(0.), cntr_Y(0.) ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF  = regn_flist[ilF] ;
    int iV0 = face_vrtx(iF,0) ;
    int iV1 = face_vrtx(iF,1) ;
    
    double sgn_F = face_regn( iF, 0 ) == iR ? +1. : -1. ;

    if ( is_curved_face(iF) ) {

      // HIGHER_ORDER approximation
      // area_X <-- +int_{[0,1]} y'(t) x(t) dt
      // area_Y <-- -int_{[0,1]} x'(t) y(t) dt
      double ax(0.), ay(0.) ; 
      double cx(0.), cy(0.) ; 
      for ( int iq=0; iq<nq; ++iq ) {
	ax += wq[iq] * face_pnts( sq[iq], iF, 0 )          * face_derv( sq[iq], iF, 1 ) ;
	ay -= wq[iq] * face_pnts( sq[iq], iF, 1 )          * face_derv( sq[iq], iF, 0 ) ;
	cx += wq[iq] * pow( face_pnts( sq[iq], iF, 0 ), 2) * face_derv( sq[iq], iF, 1 ) ;
	cy -= wq[iq] * pow( face_pnts( sq[iq], iF, 1 ), 2) * face_derv( sq[iq], iF, 0 ) ;
      }
      area_X += sgn_F * ax ;
      area_Y += sgn_F * ay ;
      cntr_X += sgn_F * cx ;
      cntr_Y += sgn_F * cy ;

    } else {

      double ome_FV = 1./2. * get_face_measure(iF) ;
      area_X += sgn_F * ome_FV * get_nor(iF,0) * ( coords_V(iV0,0) + coords_V(iV1,0) ) ;
      area_Y += sgn_F * ome_FV * get_nor(iF,1) * ( coords_V(iV0,1) + coords_V(iV1,1) ) ;
      cntr_X += sgn_F * ome_FV * get_nor(iF,0) * ( pow( coords_V(iV0,0), 2 ) + pow( coords_V(iV1,0), 2 ) + 4.*pow( coords_F(iF,0), 2 ) )/3. ;
      cntr_Y += sgn_F * ome_FV * get_nor(iF,1) * ( pow( coords_V(iV0,1), 2 ) + pow( coords_V(iV1,1), 2 ) + 4.*pow( coords_F(iF,1), 2 ) )/3. ;
  
    }
  }
  //assert( abs(area_X-area_Y)<1.e-9 ) ;
  if ( abs(area_X-area_Y)>1.e-12 ) {
    MSG("area diff --> ") ; 
    VAL(iR) ; 
    VAL(area_X) ; VAL(area_Y) ; 
    PRT( abs(abs(area_X)-abs(area_Y)) ) ;
  }

  cntr_X /= ( 2.*area_X ) ;
  cntr_Y /= ( 2.*area_Y ) ;
  
  cvd_R_coords(iR,0) = cntr_X ;
  cvd_R_coords(iR,1) = cntr_Y ;
  cvd_R_area  (iR,0) = area_X ;
}

void mesh_2Dc :: check_region_geometry() {
  MSG("begin setup_curved_regn_geom"<<endl<<flush) ;
  for ( int iR=0; iR<nR; ++iR ) {
    check_region_geometry( iR ) ;
  }
  // check
  double area_tot(0.), xbary(0.), ybary(0.) ;
  for ( int iR=0; iR<nR; ++iR ) {
    area_tot += cvd_R_area(iR,0) ;
    xbary    += cvd_R_area(iR,0) * cvd_R_coords(iR,0) ;
    ybary    += cvd_R_area(iR,0) * cvd_R_coords(iR,1) ;
  }
#if 0 // rectangular domain
  LINE(---) ;
  VAL( area_tot ) ; PRT( abs(1.-area_tot) ) ;
  VAL( xbary    ) ; PRT( abs(0.5-xbary)    ) ;
  VAL( ybary    ) ; PRT( abs(0.5-ybary)    ) ;
  LINE(---) ;
  PRT(max_diff_area) ; 
  PRT(max_diff_bary_X) ; 
  PRT(max_diff_bary_Y) ; 
  LINE(---) ;
#else // circular domain
  LINE(---) ;
  double pi = acos(-1.) ;
  VAL( area_tot ) ; PRT( abs(pi-abs(area_tot)) ) ;
  VAL( xbary    ) ; PRT( ybary ) ;
  LINE(---) ;
  PRT(max_diff_area) ; 
  PRT(max_diff_bary_X) ; 
  PRT(max_diff_bary_Y) ; 
  LINE(---) ;
#endif
  MSG("end-->setup_curved_regn_geom"<<endl<<flush) ;
}

void mesh_2Dc :: check_region_geometry( int iR ) { 
  // this routine works nicely if R is star-shaped with
  // respect to the point (xE,yE) given by the arithmetic mean.
  assert(DIM==2) ;

  // get region vertices
  vector<int> vlist ;
  get_regn_vrtx( iR, vlist ) ;
  const int nRV = vlist.size() ;
    
  double xV[nRV], yV[nRV] ;
  double xE(0.), yE(0.) ;
  for ( int ilV=0; ilV<nRV; ++ilV ) {
    int iV = vlist[ilV] ;
    xV[ilV] = coords_V(iV,0) ;
    yV[ilV] = coords_V(iV,1) ;
    xE += xV[ilV] ;
    yE += yV[ilV] ;
  }
  xE /= double(nRV) ; 
  yE /= double(nRV) ; 

  // get region vertices
  vector<int> flist ;
  get_regn_face( iR, flist ) ;
  const int nRF = flist.size() ;
  assert( nRF==nRV ) ;
  
  // loop on the edges of the element
  double xR(0.), yR(0.), aR(0.) ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = flist[ilF] ;
    if ( is_curved_face(iF) ) {

      int max_np = 10 ; // 10 is a safe factor
      for ( int ip=0; ip<max_np; ++ip ) {
	
	double t0 = double(ip)  /double(max_np) ;
	double t1 = double(ip+1)/double(max_np) ;
	
	double xq[DIM+1] = { xE, face_pnts(t0,iF,0), face_pnts(t1,iF,0) } ;
	double yq[DIM+1] = { yE, face_pnts(t0,iF,1), face_pnts(t1,iF,1) } ;
	double aR_Ti = abs( (xq[0]-xq[2])*(yq[0]-yq[1])-(yq[0]-yq[2])*(xq[0]-xq[1]) )/2. ;
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

    } else {

      // quadrature stuff on the sub-triangle Ti (ilV,ilV+1,iR)
      // (vrtx quad rule exact for linear functions) 
      int ilV = ilF ;
      int il0 = ilV ;
      int il1 = (ilV+1)%nRV ;
      double xq[DIM+1] = { xE, xV[il0], xV[il1] } ;
      double yq[DIM+1] = { yE, yV[il0], yV[il1] } ;
      double aR_Ti = abs( (xq[0]-xq[2])*(yq[0]-yq[1])-(yq[0]-yq[2])*(xq[0]-xq[1]) )/2. ;
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

    } // end of    if ( is_curved_face(iF) ) {... } else {...
  } // end of   for ( int ilF=0; ilF<nRF; ++ilF ) {...

  // -----------------------------------------------
  max_diff_area   = max( max_diff_area,   abs(aR   -cvd_R_area  (iR,0) ) ) ;
  max_diff_bary_X = max( max_diff_bary_X, abs(xR/aR-cvd_R_coords(iR,0) ) ) ;
  max_diff_bary_Y = max( max_diff_bary_Y, abs(yR/aR-cvd_R_coords(iR,1) ) ) ;
  // -----------------------------------------------
}

void mesh_2Dc :: print_region_data() {
  for ( int iR=0; iR<nR; ++iR ) {
    print_region_data( iR ) ;
  }
}
void mesh_2Dc :: print_region_data( int iR ) {
  const string face_type_str[2] = { "RADIAL FACE",
				    "AZIMUTH FACE" } ;
  LINE(--) ;
  MSG("print region data --> ") ; PRT(iR) ; 
  vector<int> cvd_regn_flist ;
  get_regn_cvd_face( iR, cvd_regn_flist ) ;
  int ncF = cvd_regn_flist.size() ;
  for ( int ilF=0; ilF<ncF; ++ilF ) {
    int iF = cvd_regn_flist[ilF] ;
    cout << "ilF = " << ilF << ", iF = " << iF << endl ;
    cout << "    gl    " << face_type_str[ cvd_flist[iF].get_grid_line() ] << endl ;
    cout << "    gv    " << cvd_flist[iF].get_grid_val() << endl ;
    cout << "    s0,s1 " << cvd_flist[iF].get_s0()  << ", " << cvd_flist[iF].get_s1() << endl ;
    double x0 = face_pnts( 0., iF, 0 ) ;
    double y0 = face_pnts( 0., iF, 1 ) ;
    double x1 = face_pnts( 1., iF, 0 ) ;
    double y1 = face_pnts( 1., iF, 1 ) ;
    cout << "    (" << x0 << "," << y0 << ") --> (" << x1 << "," << y1 << ")" << endl ;
    double xm_s = face_derv( 0.5, iF, 0 ) ;
    double ym_s = face_derv( 0.5, iF, 1 ) ;
    cout << "    der (" << xm_s << "," << ym_s << ") " << endl ;
    double eps = 1.e-8 ;
    double dder = 
      ( face_pnts(0.5+eps,iF,0) - face_pnts(0.5-eps,iF,0) )/( 2.*eps ) ;
    cout << "    dder = " << dder << endl ;
  }
}

#endif // end of _MESH2DC_HH
