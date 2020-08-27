#ifndef _MESH_2D_BUILTIN_HH
#define _MESH_2D_BUILTIN_HH

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// COORDINATE MAPS
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class KershawMap {
private:
  // input members
  int nx    ; // number of partitions along X
  int ny    ; // number of partitions along Y
  double wi ; // mapping parameter along X
  double wj ; // mapping parameter along Y

  // local vars
  double dx, dy0, dy1 ;
  int    kx, ky ;

  // mapping aux functions
  int rflag( int i ) const {
    int ival(0) ;
    if      (       0<=i && i<kx      ) { ival = 0 ; }
    else if (      kx<=i && i<2*kx    ) { ival = 1 ; }
    else if (    2*kx<=i && i<nx-2*kx ) { ival = 2 ; }
    else if ( nx-2*kx<=i && i<nx-kx   ) { ival = 3 ; }
    else if ( nx  -kx<=i && i<=nx     ) { ival = 4 ; }
    else { assert(false) ; }
    return ival ;
  }
  double y_left( int j ) const {
    double retval(0.) ;
    if ( 0<=j && j<ky ) { retval = dy0*double(j) ; } 
    else                { retval = wj+dy1*double(j-ky) ; } 
    return retval ;
  }
  double y_right( int j ) const {
    double retval(0.) ;
    if ( 0<=j && j<ky ) { retval = dy1*double(j) ; } 
    else                { retval = (1.-wj)+dy0*double(j-ky) ; } 
    return retval ;
  }
  // straight lines for Kershaw maps
  double r0( int i, int j ) const {
    assert( 0<=i && i<kx ) ;
    return y_left (j) ; 
  }
  double r1( int i, int j ) const { 
    assert( kx<=i && i<2*kx ) ;
    return y_left(j)+( y_right(j)-y_left(j) ) * double(i-kx)/double(kx) ; 
  }
  double r2( int i, int j ) const {
    assert( 2*kx<=i && i<nx-2*kx ) ;
    return y_right(j)+( y_left(j)-y_right(j) ) * double(i-2*kx)/double(nx-4*kx) ; 
  }
  double r3( int i, int j ) const {
    assert( nx-2*kx<=i && i<nx-kx ) ;
    return y_left(j)+( y_right(j)-y_left(j) ) * double(i-nx+2*kx)/double(kx) ; 
  }
  double r4( int i, int j ) const { 
    assert( nx-kx<=i && i<=nx ) ;
    return y_right(j) ; 
  }
  
public:
  KershawMap( int _nx=2, int _ny=2, double _wi=0.25, double _wj=0.5 ) : 
    nx(_nx), ny(_ny), wi(_wi), wj(_wj) {

    assert( 0.<wi && wi<=0.25 ) ;
    assert( 0.<wj && wj<=0.5  ) ;
    
    assert( nx>0 && nx%2==0 ) ;
    assert( ny>0 && ny%2==0 ) ;
    
    kx = int( wi*double(nx) ) ;
    ky = ny/2 ;
    
    dx  = 1./double(nx) ;

    dy0 = wj/double(ky) ;
    dy1 = (1.-wj)/double(ky) ;
  }
  ~KershawMap() {}

  // ---------------
  inline double xij( int i, int j ) const {
    return dx * double(i) ;
  }
  inline double yij( int i, int j ) const { 
    double retval(0.) ;
    switch( rflag(i) ) {
    case 0: retval = r0(i,j) ; break ;
    case 1: retval = r1(i,j) ; break ;
    case 2: retval = r2(i,j) ; break ;
    case 3: retval = r3(i,j) ; break ;
    case 4: retval = r4(i,j) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
} ;

class CoordinateMap {
private:
  // input members
  int nx    ; // number of partitions along X
  int ny    ; // number of partitions along Y
  int mflag ; // mesh flag to select the mapping function
  double wi ; // mapping parameter along X
  double wj ; // mapping parameter along Y

  // local vars
  double pi, dx, dy ;

  // pointer for Kershaw Mapping
  KershawMap * p_kmap ;

  // aux methods
  inline double xy_rand( double a=-1., double b=1. ) const {
    //return a + (b-a) * double(random())/double(LONG_MAX) ;
    return a + (b-a) * double(random())/double(RAND_MAX) ;
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
  inline bool is_equal( double x, double y ) const {
    return abs(x-y)<1.e-13 ;
  }

#if 0
  // two rectangular regions inside, [0.2,0.25]x[0,0.6] and [0.7,0.75]x[0.4,1] and
  inline bool is_bnd_vertex1( int i, int j ) const {
    return is_bnd_vertex(i,j) ||
      ( ( is_equal(dx*i,0.2) || is_equal(dx*i,0.25) ) && ( dy*j<0.6 || is_equal(dy*j,0.6) ) ) ||
      ( ( is_equal(dx*i,0.7) || is_equal(dx*i,0.75) ) && ( dy*j>0.4 || is_equal(dy*j,0.4) ) ) ;
  }

  inline bool is_bnd_vertex1( double x, double y ) const {
    int i = int(x/dx) ;
    int j = int(y/dy) ;
    return is_bnd_vertex(i,j) ||
      ( ( is_equal(dx*i,0.2) || is_equal(dx*i,0.25) ) && ( dy*j<0.6 || is_equal(dy*j,0.6) ) ) ||
      ( ( is_equal(dx*i,0.7) || is_equal(dx*i,0.75) ) && ( dy*j>0.4 || is_equal(dy*j,0.4) ) ) ;
  }
#endif

  // vertices on line i=nx/2 are considered as internal boundary
  // (used to randomize a mesh with a fixed internal line)
  inline bool is_bnd_vertex1( int i, int j ) const {
    return is_bnd_vertex(i,j) || is_equal(dx*i,0.5) ;
  }
  inline bool is_bnd_vertex1( double x, double y ) const {
    int i = int(x/dx) ;
    int j = int(y/dy) ;
    return is_bnd_vertex(i,j) || is_equal(x,0.5) ;
  }

  // map-0, identity
  // ---------------
  inline double s0( int i, int j )       const { return dx * double(i) ; }
  inline double t0( int i, int j )       const { return dy * double(j) ; }
  inline double s0( double x, double y ) const { return x ; }
  inline double t0( double x, double y ) const { return y ; }

#if 0 // modified for the paper with Daniil on KKT
  // map-1 (sin(..)sin(..)-Los Alamos map)
  // ----------------------------------
  inline double s1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    double cij = 50.*tij*(1.-tij)*(3./8.-tij)*(5./8.-tij) ;
    //double cij = (tij-3./8.)*(5./8.-tij) < 0. ? 
    //  50.*tij*(1.-tij)*(3./8.-tij)*(5./8.-tij) : 0. ;
    return sij + wi *cij*sin(2.*pi*sij) ; 
  }
  inline double t1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    double cij = 50.*tij*(1.-tij)*(3./8.-tij)*(5./8.-tij) ;
    //double cij = (tij-3./8.)*(5./8.-tij) < 0. ? 
    //  50.*tij*(1.-tij)*(3./8.-tij)*(5./8.-tij) : 0. ;
    return tij + wj * cij*sin(2.*pi*sij) ; 
  }
  inline double s1( double x, double y ) const { 
    return x + wi * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1( double x, double y ) const { 
    return y + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
#endif

#if 1 // ORIGINAL
  // map-1 (sin(..)sin(..)-Los Alamos map)
  // ----------------------------------
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
  inline double s1( double x, double y ) const { 
    return x + wi * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1( double x, double y ) const { 
    return y + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
#endif

#if 0 // top-right quadrant, boundary cells are kept fixed
  inline double s1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    bool no_action = (i==1 || i==nx-1 || j==1 || j==ny-1) ;
    return no_action ? sij :
      sij + wi * sin(pi*sij)*sin(pi*tij) ; 
  }
  inline double t1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    bool no_action = (i==1 || i==nx-1 || j==1 || j==ny-1) ;
    return no_action ? tij :
      tij + wj * sin(pi*sij)*sin(pi*tij) ; 
  }
  inline double s1( double x, double y ) const { 
    return x + wi * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1( double x, double y ) const { 
    return y + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
#endif

  // !!!!
#if 0 // top-right quadrant, boundary cells are moving
  // map-1x (sin(..)sin(..)-Los Alamos map)
  // ----------------------------------
  inline double s1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    double pij = sin(pi*sij)*sin(pi*tij) ; 
    double wij = wi * ( 1.-exp( pij ) ) ;
    return sij + wij ;
  }
  inline double t1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    double pij = sin(pi*sij)*sin(pi*tij) ; 
    double wij = wj * ( 1.-exp( pij ) ) ;
    return tij + wij ;
  }
  inline double s1( double x, double y ) const { 
    return x + wi * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1( double x, double y ) const { 
    return y + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
#endif


#if 0 // different mapping, equivalent to top-right quadrant
  // map-1b (sin(..)sin(..)-Los Alamos map)
  // ----------------------------------
  inline double s1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return sij + wi * 4096./9. * sij*(1.-sij)*(1./2.-sij)*tij*(1.-tij)*(1./2.-tij) ; 
  }
  inline double t1( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return tij + wj * 4096./9. * sij*(1.-sij)*(1./2.-sij)*tij*(1.-tij)*(1./2.-tij) ; 
  }
  inline double s1( double x, double y ) const { 
    return x + wi * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
  inline double t1( double x, double y ) const { 
    return y + wj * sin(2.*pi*x)*sin(2.*pi*y) ; 
  }
#endif

  // map-2 (randomized vertices)
  // ---------------------------
  inline double s2( int i, int j ) const { 
    double sij = s0(i,j) ; 
    //a + (b-a) * double(random())/double(RAND_MAX) ;
    return is_bnd_vertex(i,j) ? sij : sij + wi * xy_rand( -dx/2., dx/2. ) ;
  }
  inline double t2( int i, int j ) const { 
    double tij = t0(i,j) ; 
    return is_bnd_vertex(i,j) ? tij : tij + wj * xy_rand( -dy/2., dy/2. ) ;
  }
  inline double s2( double x, double y ) const { 
    return x + wi * xy_rand( -dx/2., dx/2. ) ;
  }
  inline double t2( double x, double y ) const { 
    return y + wj * xy_rand( -dy/2., dy/2. ) ;
  }

  // map-21 (randomized vertices with two regions inside)
  // ---------------------------
  inline double s21( int i, int j ) const { 
    double sij = s0(i,j) ; 
    return is_bnd_vertex1(i,j) ? sij : sij + wi * xy_rand( -dx/2., dx/2. ) ;
  }
  inline double t21( int i, int j ) const { 
    double tij = t0(i,j) ; 
    return is_bnd_vertex1(i,j) ? tij : tij + wj * xy_rand( -dy/2., dy/2. ) ;
  }
  inline double s21( double x, double y ) const { 
    //assert(false) ; // do not use this method
    return is_bnd_vertex1(x,y) ? x : x + wi * xy_rand( -dx/2., dx/2. ) ;
  }
  inline double t21( double x, double y ) const { 
    //assert(false) ; // do not use this method
    return is_bnd_vertex1(x,y) ? y : y + wj * xy_rand( -dy/2., dy/2. ) ;
  }

  // map-3 (distorted quad)
  // ----------------------
  inline double s3( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return sij + wi * tij ;
  }
  inline double t3( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ; 
    return wj *sij + tij ;
  }
  inline double s3( double x, double y ) const { 
    return x + wi * y ;
  }
  inline double t3( double x, double y ) const { 
    return wj * y + x ;
  }

  // map-4 (Kershaw's mesh)
  // ----------------------
  inline double s4( int i, int j ) const { return p_kmap->xij(i,j) ; }
  inline double t4( int i, int j ) const { return p_kmap->yij(i,j) ; }

  // map-5 (distorted quad, second implementation, constant aspect ratio)
  // ----------------------
  // wi == 0 means "rectangle"
  // wi is the angle between axis Y and the oblique edge 
  // of the parallelogram (also wi = 90-theta)
  // wj is used as aspect ratio, wi>1 means Y-side > X-side
  inline double s5( int i, int j ) const {
    assert( wj>0. ) ;
    double sij = s0(i,j) ; 
    double tij = wj * t0(i,j) ; 
    return sij + cos(pi/2-pi*wi/180.) * tij ;
  }
  inline double t5( int i, int j ) const { 
    assert( wj>0. ) ;
    double sij = s0(i,j) ; 
    double tij = wj * t0(i,j) ; 
    return 0.  + sin(pi/2-pi*wi/180.) * tij ;
  }
  inline double s5( double x, double y ) const { 
    assert( wj>0. ) ;
    return x  + cos(pi/2-pi*wi/180.) * y * wj ;
  }
  inline double t5( double x, double y ) const { 
    assert( wj>0. ) ;
    return 0. + sin(pi/2-pi*wi/180.) * y * wj ;
  }

  // map-6 distorted mesh (see Ciarlet Theorem for FEM)
  // --------------------------------------------------
  inline double s6( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ;
    assert( sij>=0. && wi>0. ) ; 
    return pow( abs(sij), wi ) ;
  }
  inline double t6( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ;
    assert( tij>=0. && wj>0. ) ;  
    return pow( abs(tij), wj ) ;
  }
  inline double s6( double x, double y ) const { 
    assert( x>=0. && wi>0. ) ;  
    return pow( x, wi ) ;
  }
  inline double t6( double x, double y ) const { 
    assert( y>=0. && wj>0. ) ;  
    return pow( y, wj ) ;
  }

  // map-8 (slope quad)
  // ----------------------
  // w0 = wi, w1 = wj
  inline double s8( int i, int j ) const { 
    double sij = s0(i,j) ; 
    return sij ;
  }
  inline double t8( int i, int j ) const { 
    double sij = s0(i,j) ; 
    double tij = t0(i,j) ;
    double w0(wi), w1(wj) ;
    double t0 = w0 + double(j)/double(ny)*(1.-w0) ;
    double t1 = double(j)/double(ny)*(1.-w1) ;
    return t0 + (t1-t0) * sij ;
  }
  inline double s8( double x, double y ) const { 
    return x ;
  }
  inline double t8( double x, double y ) const {
    double w0(wi), w1(wj) ;
    double t0 = w0 + y*(1.-w0) ;
    double t1 = y*(1.-w1) ;
    return t0 + (t1-t0) * x ;
  }

  // map-9 (slope quad)
  // ----------------------
  // w0 = wi, w1 = wj
  inline double s9( int i, int j ) const { 
    bool ok_midcell = 
      ( i==nx/2   && j==ny/2   ) || 
      ( i==nx/2+1 && j==ny/2   ) ||
      ( i==nx/2+1 && j==ny/2+1 ) ||
      ( i==nx/2   && j==ny/2+1 ) ;
    bool ok_column = i==(3*nx)/4 || i==(3*nx)/4+1 ;
    double sij = ok_midcell || ok_column ?
      s0(i,j) : s2(i,j) ; 
    return sij ;
  }
  inline double t9( int i, int j ) const { 
    bool ok_midcell = 
      ( i==nx/2   && j==ny/2   ) || 
      ( i==nx/2+1 && j==ny/2   ) ||
      ( i==nx/2+1 && j==ny/2+1 ) ||
      ( i==nx/2   && j==ny/2+1 ) ;
    bool ok_column = i==(3*nx)/4 || i==(3*nx)/4+1 ;
    double tij = ok_midcell || ok_column ?
      t0(i,j) : t2(i,j) ; 
    return tij ;
  }
  inline double s9( double x, double y ) const { 
    return s9( int(x/dx), int(y/dy) ) ;
  }
  inline double t9( double x, double y ) const {
    return t9( int(x/dx), int(y/dy) ) ;
  }

  
public:
  CoordinateMap( int _nx=1, int _ny=1, int _mflag=0, double _wi=0., double _wj=0. ) : 
    nx(_nx), ny(_ny), mflag(_mflag), wi(_wi), wj(_wj) {
    assert( nx>0 && ny>0 ) ;
    dx = 1./double(nx) ;
    dy = 1./double(ny) ;
    pi = acos(-1.) ;
    srandom(2) ;
    if ( mflag==4 ) {
      p_kmap = new KershawMap(nx,ny,wi,wj) ;
    }
  }
  ~CoordinateMap() {}
  // set up methods
  inline void set_wij( double _wi, double _wj ) {
    wi=_wi ;
    wj=_wj ;
  }
  // return vrtx index
  inline int vrtx_index( int i, int j ) const {
    return i + (nx+1) * j ;
  }
  // map driver methods
  inline double sxy( int i, int j ) const {
    double retval(0.) ;
    switch( mflag ) {
    case  0: retval=s0(i,j) ; break ;
    case  1: retval=s1(i,j) ; break ;
    case  2: retval=s2(i,j) ;  break ;
    case  3: retval=s3(i,j) ; break ;
    case  4: retval=s4(i,j) ; break ;
    case  5: retval=s5(i,j) ; break ;
    case  6: retval=s6(i,j) ; break ;
    case  7: retval=s21(i,j) ; break ;
    case  8: retval=s8(i,j) ; break ;
    case  9: retval=s9(i,j) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy( int i, int j ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case  0: retval=t0(i,j) ; break ;
    case  1: retval=t1(i,j) ; break ;
    case  2: retval=t2(i,j) ; break ;
    case  3: retval=t3(i,j) ; break ;
    case  4: retval=t4(i,j) ; break ;
    case  5: retval=t5(i,j) ; break ;
    case  6: retval=t6(i,j) ; break ;
    case  7: retval=t21(i,j) ; break ;
    case  8: retval=t8(i,j) ; break ;
    case  9: retval=t9(i,j) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double sxy( double x, double y ) const {
    double retval(0.) ;
    switch( mflag ) {
    case  0: retval=s0(x,y) ; break ;
    case  1: retval=s1(x,y) ; break ;
    case  2: retval=s2(x,y) ; break ;
    case  3: retval=s3(x,y) ; break ;
    case  4: MSG("no Kershaw's mesh") ; assert(false)  ; break ;
    case  5: retval=s5(x,y) ; break ;
    case  6: retval=s6(x,y) ; break ;
    case  7: retval=s21(x,y) ; break ;
    case  8: retval=s8(x,y) ; break ;
    case  9: retval=s9(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy( double x, double y ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case  0: retval=t0(x,y) ; break ;
    case  1: retval=t1(x,y) ; break ;
    case  2: retval=t2(x,y) ; break ;
    case  3: retval=t3(x,y) ; break ;
    case  4: MSG("no Kershaw's mesh") ; assert(false)  ; break ;
    case  5: retval=t5(x,y) ; break ;
    case  6: retval=t6(x,y) ; break ;
    case  7: retval=t21(x,y) ; break ;
    case  8: retval=t8(x,y) ; break ;
    case  9: retval=t9(x,y) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }

  void bil_map_xy( double xh, double yh, vector<double> & xvrt, vector<double> & yvrt, double & x, double & y ) {
    assert( xvrt.size()==4  && yvrt.size()==4  ) ;
    assert( 0<=xh && xh<=1. && 0<=yh && yh<=1. ) ;
    vector<double> phi_vec(4) ;
    phi_vec[0] = (1.-xh)*(1.-yh) ;
    phi_vec[1] = xh*(1.-yh) ;
    phi_vec[2] = xh*yh ;
    phi_vec[3] = (1.-xh)*yh ;
    x = y = 0. ;
    for ( int il=0; il<4; ++il ) {
      x += phi_vec[il] * xvrt[il] ;
      y += phi_vec[il] * yvrt[il] ;
    }
  }
} ;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// QUADRILATERAL-BASED MESH (standard flag: 1xx)
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void build_quad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
		      double wx=0., double wy=0. ) {
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
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      int iV = cmap.vrtx_index(i,j) ;
      xV[iV] = cmap.sxy(i,j) ;
      yV[iV] = cmap.txy(i,j) ;
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
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
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
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}

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
#if 0
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
      if ( mesh_flag/10==1 || (mesh_flag/10==2 && i%2==0) ) { // 11?, by translation
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
      } else { // 12? by reflection
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
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("JQuad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
#endif
//#include "jquad.hh"
//-----------------------------------------------------------------------------------------
// modified quad mesh, corner quadrilaterals (two edges on the boundary) are split into 
// triangular cells
//
//     |-----------|
//     |\ |  |  | /|
//     |-----------|
//     |  |  |  |  |
//     |-----------|
//     |  |  |  |  |
//     |-----------|
//     |/ |  |  | \|
//     |-----------|
//
//     (example of 2x2 mesh)
//-----------------------------------------------------------------------------------------
void build_mquad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
		      double wx=0., double wy=0. ) {
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
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      int iV = cmap.vrtx_index(i,j) ;
      xV[iV] = cmap.sxy(i,j) ;
      yV[iV] = cmap.txy(i,j) ;
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
      bool corner_cell = 
	bool(i==0    && j==0) || bool(i==0    && j==ny-1) || 
	bool(i==nx-1 && j==0) || bool(i==nx-1 && j==ny-1) ;
      if ( !corner_cell ) {
	regn_vlist.push_back( 4 ) ;
	regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
	regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
	regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
	regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
	fR.push_back( UNSET ) ;
      }
    }
  }
  { // i==0 && j==0
    int i(0), j(0) ;
    // bottom cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    fR.push_back( UNSET ) ;
    // top cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    fR.push_back( UNSET ) ;
  }
  { // i==0 && j==ny-1
    int i(0), j(ny-1) ;
    // bottom cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    fR.push_back( UNSET ) ;
    // top cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    fR.push_back( UNSET ) ;
  }
  { // i==nx-1 && j==0 
    int i(nx-1), j(0) ;
    // bottom cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    fR.push_back( UNSET ) ;
    // top cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    fR.push_back( UNSET ) ;
  }
  { // i==nx-1 && j==ny-1
    int i(nx-1), j(ny-1) ;
    // bottom cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,  j) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    fR.push_back( UNSET ) ;    
    // top cell
    regn_vlist.push_back( 3 ) ;
    regn_vlist.push_back( cmap.vrtx_index(i+1,j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,  j+1) ) ;
    regn_vlist.push_back( cmap.vrtx_index(i,    j) ) ;
    fR.push_back( UNSET ) ;    
  }
  
  int nR = nx*ny+4 ;
  regn_vlist.push_back( nR ) ;
  VAL(nR) ; VAL(regn_vlist.size()) ; PRT(fR.size()) ;
  assert( (regn_vlist.size()-8*4)/5+8==fR.size() ) ;
  assert( nR==(regn_vlist.size()-8*4)/5+8 ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
//#include "mesh2D_quadhole.hh"
#include "mesh2D_circleframe.hh"
//#include "mesh2D_two_skewed_quads.hh"
#if 1
void build_new_quad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
			  double wx=0., double wy=0. ) {
  switch( mesh_flag/10 ) {
  case 0: build_quad_mesh ( mesh, nx, ny, mesh_flag, wx, wy ) ; break ;
  case 1: build_cquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ;
  //case 1: build_jquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ; // jerome meshes
  //case 2: build_jquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ; // jerome meshes
  //case 3: build_jquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ; // jerome meshes with double reflection
  //case 4: build_hquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ; // central hole (monotonicity for MFD)
  //case 5: build_mquad_mesh( mesh, nx, ny, mesh_flag, wx, wy ) ; break ; // central hole (monotonicity for MFD)
  default: assert(false) ;
  }
}
#endif
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// REGULAR-SHAPED MESH OF TRIANGLES (standard flag: 2xx)
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*
 top-right               top-left                          top          
 (i,j+1)  (i+1,j+1)      (i,j+1)  (i+1,j+1)           (i,j+1)  (i+1,j+1)
  |--------|              |--------|                   |--------|
  |\       |              |       /|                   |\      /|
  | \      |              |      / |                   | \    / |
  |  \     |              |     /  |                   |  \  /  |
  |   \    |              |    /   |             left  |   \/   |  right
  |    \   |              |   /    |                   |   /\   | 
  |     \  |              |  /     |                   |  /  \  |
  |      \ |              | /      |                   | /    \ |
  |       \|              |/       |                   |/      \|  
  |--------|              |--------|                   |--------|
 (i,j)    (i+1,j)        (i,j)    (i+1,j)             (i,j)    (i+1,j)
 bottom-left             bottom-right                    bottom
*/ 
int vrtx_indx( int i, int j, int nx ) {
  return i + j * (nx+1) ;
}
int cell_indx( int i, int j, int nx, int nV ) {
  return nV + i + j * nx ;
}
//-------------------------------------------------------------------------------
// diagonal from top-left to bottom-right corner
//-------------------------------------------------------------------------------
void bottom_left_triangle ( int i, int j, int nx, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( vrtx_indx(i,  j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j+1, nx ) ) ;
}
void top_right_triangle   (  int i, int j, int nx, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( vrtx_indx(i,  j+1, nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j+1, nx ) ) ;
}
//-------------------------------------------------------------------------------
// diagonal from top-right to bottom-left corner
//-------------------------------------------------------------------------------
void bottom_right_triangle( int i, int j, int nx, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( vrtx_indx(i,  j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j+1, nx ) ) ;
}
void top_left_triangle    ( int i, int j, int nx, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( vrtx_indx(i,  j,   nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j+1, nx ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j+1, nx ) ) ;
}
//-------------------------------------------------------------------------------
// criss-cross diagonals
//-------------------------------------------------------------------------------
void bottom_triangle ( int i, int j, int nx, int nV, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( cell_indx(i,  j,   nx, nV ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j,   nx     ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j,   nx     ) ) ;
}
void top_triangle    ( int i, int j, int nx, int nV, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( cell_indx(i,  j,   nx, nV ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j+1, nx     ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j+1, nx     ) ) ;
}
void right_triangle  ( int i, int j, int nx, int nV, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( cell_indx(i,  j,   nx, nV ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j,   nx     ) ) ;
  regn_vlist.push_back( vrtx_indx(i+1,j+1, nx     ) ) ;
}
void left_triangle   ( int i, int j, int nx, int nV, vector<int> & regn_vlist ) {
  regn_vlist.push_back(3) ;
  regn_vlist.push_back( cell_indx(i,  j,   nx, nV ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j+1, nx     ) ) ;
  regn_vlist.push_back( vrtx_indx(i,  j,   nx     ) ) ;

}
//-------------------------------------------------------------------------------
void build_tria_mesh( mesh_2Dv & mesh, 
		      int    nx=1,  int    ny=1,   int mesh_flag=0,
		      double wx=0., double wy=0. ) {
  assert( nx>0 && ny>0 ) ;
  // ----------------------------------------------------------------------------
  // set vertex coordinates
  // ----------------------------------------------------------------------------
  int cmap_flag = mesh_flag % 10 ;
  int topl_flag = ( mesh_flag / 10 ) % 10 ;
  VAL(mesh_flag) ;
  VAL(cmap_flag) ;  
  PRT(topl_flag) ;
  CoordinateMap cmap(nx,ny,cmap_flag,wx,wy) ;
  vector<double> xV, yV ;
  for ( int j=0; j<=ny; ++j ) {
    for ( int i=0; i<=nx; ++i ) {
      double xij = cmap.sxy(i,j) ;
      double yij = cmap.txy(i,j) ;
      xV.push_back(xij) ;
      yV.push_back(yij) ;
    }
  }
  if ( topl_flag==2 ) {
    for ( int j=0; j<ny; ++j ) {
      for ( int i=0; i<nx; ++i ) {
	int iV0 = vrtx_indx( i,   j,   nx ) ;
	int iV1 = vrtx_indx( i+1, j,   nx ) ;
	int iV2 = vrtx_indx( i+1, j+1, nx ) ;
	int iV3 = vrtx_indx( i,   j+1, nx ) ;
	double xij = ( xV[iV0]+xV[iV1]+xV[iV2]+xV[iV3] )/4. ;
	double yij = ( yV[iV0]+yV[iV1]+yV[iV2]+yV[iV3] )/4. ;
	xV.push_back(xij) ;
	yV.push_back(yij) ;
      }
    }
    int n_vrtx = xV.size() ;
    PRT(n_vrtx) ;
    assert( xV.size()==yV.size() ) ;
  }
  // ----------------------------------------------------------------------------
  // set region connectivity
  // ----------------------------------------------------------------------------
  int nR = 0 ;
  int nV = (nx+1)*(ny+1) ;
  vector<int> regn_vlist ;
  switch( topl_flag ) {
  case 0:
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j ) {
	bottom_left_triangle(i,j,nx,regn_vlist) ;
	top_right_triangle  (i,j,nx,regn_vlist) ;
      }
    }
    nR = 2*nx*ny ;
    break ;
  case 1:
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j ) {
	bottom_right_triangle(i,j,nx,regn_vlist) ;
	top_left_triangle    (i,j,nx,regn_vlist) ;
      }
    }
    nR = 2*nx*ny ;
    break ;
  case 2:
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j ) {
	bottom_triangle(i,j,nx,nV,regn_vlist) ;
	top_triangle   (i,j,nx,nV,regn_vlist) ;
	left_triangle  (i,j,nx,nV,regn_vlist) ;
	right_triangle (i,j,nx,nV,regn_vlist) ;
      }
    }
    nR  = 4*nx*ny ;
    nV += nx*ny ;
    break ;
  default :
    assert(0) ;
  }
  regn_vlist.push_back(nR) ;
  assert(nR>0) ;
  //-----------------------------------------------------------------------------------------
  vector<int> fV(nV), fR(nR) ;
  for ( int iV=0; iV<nV; ++iV ) { fV[iV]=UNSET ; }
  for ( int iR=0; iR<nR; ++iR ) { fR[iR]=UNSET ; }
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR REGULAR TRIANGULAR MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 2xx)
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",200+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("Tria-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// L-SHAPED (NON-CONVEX) CELL CONSTRUCTION (standard flag: 3xx)
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*

    4        5        6        7
    o--------o--------o--------o
    |        |                 |
    |        |      cell 1     |
    |        |                 |
    |        |                 |
    |        o--------o        |
    |        8       9|        |
    |                 |        |
    |    cell 0       |        |
    |                 |        |
    o--------o--------o--------o
    0        1        2        3

    (local numbering)

*/
//-----------------------------------------------------------------------------------------
class LShapeMesh {
private:
  int nx, ny, mesh_flag ;
  double wx, wy ;
  int ** vmap ;
  static const int vind[2][7] ;
  
  void allocate_vmap() ;
  int vrtx_index( int ir, int il, int I, int J ) ;

public:
  LShapeMesh( int _nx=1, int _ny=1, int _mesh_flag=0, double _wx=0., double _wy=0. ) : 
    nx(_nx), ny(_ny), mesh_flag(_mesh_flag), wx(_wx), wy(_wy) {}
  ~LShapeMesh() {}

  void build_region_topology( vector<int> & regn_vlist, vector<int> & fR ) ;
  void build_vertex_coords  ( vector<double> & xV, vector<double> & yV, vector<int> & fV ) ;
} ;

const int LShapeMesh::vind[2][7] = 
  { 
    { 0, 1, 2, 9, 8, 5, 4 },
    { 2, 3, 7, 6, 5, 8, 9 }
  } ;

void LShapeMesh::allocate_vmap() {
  const int nvx = 3*nx+1 ;
  const int nvy = 2*ny+1 ;
  try {
    vmap = new int*[nvx] ;
    for ( int i=0; i<nvx; ++i ) {
      vmap[i] = new int[nvy] ;
    }
  } catch(...) {
    cerr << "LshapeMap: error in memory allocation of vmap" << endl ;
    exit(0) ;
  }
  // if OK, then reset to UNSET
  for ( int j=0; j<nvy; ++j ) {
    for ( int i=0; i<nvx; ++i ) {
      vmap[i][j] = UNSET ;
    }
  }
}

int LShapeMesh::vrtx_index( int ir, int il, int I, int J ) {
  assert( ir==0 || ir==1 ) ;
  assert( 0<=il && il<7  ) ;
  assert( 0<=I  && I<nx  ) ;
  assert( 0<=J  && J<ny  ) ;
  int k = vind[ir][il] ;
  int i = k<8 ? 3*I + k % 4   : 3*I + k-7 ;
  int j = k<8 ? 2*J + 2*(k/4) : 2*J + 1   ;
  assert( i<=3*nx+1 && j<=2*ny+1 ) ;
  return vmap[i][j] ;
}

void LShapeMesh::build_vertex_coords( vector<double> & xV, vector<double> & yV, vector<int> & fV ) {
  assert( xV.size()==0 && yV.size()==0 && fV.size()==0 ) ;
  
  // allocate and reset vmap before building vertices
  allocate_vmap() ;
  
  // set total number of vertices
  const int nvx = 3*nx+1 ;
  const int nvy = 2*ny+1 ;

  CoordinateMap cmap(nvx-1,nvy-1,mesh_flag,wx,wy) ;

  // build "external" L-shape vertices
  int iV = 0;
  for ( int j=0; j<nvy; j+=2 ) {
    for ( int i=0; i<nvx; ++i ) {
      xV.push_back( cmap.sxy(i,j) ) ;
      yV.push_back( cmap.txy(i,j) ) ;
      fV.push_back( UNSET ) ;
      vmap[i][j] = iV++ ;
    }
  }
  // build "internal" L-shape vertices
  for ( int j=1; j<nvy; j+=2 ) {
    for ( int k=0; k<nx; ++k ) {
      { // first internal L-shape vertex
	int i  = 3*k+1 ;
	xV.push_back( cmap.sxy(i,j) ) ;
	yV.push_back( cmap.txy(i,j) ) ;
	fV.push_back( UNSET ) ;
	vmap[i][j] = iV++ ;
      }
      { // second internal L-shape vertex
	int i  = 3*k+2 ;
	xV.push_back( cmap.sxy(i,j) ) ;
	yV.push_back( cmap.txy(i,j) ) ;
	fV.push_back( UNSET ) ;
	vmap[i][j] = iV++ ;
      }
    }
  }
}

void LShapeMesh::build_region_topology( vector<int> & regn_vlist, vector<int> & fR ) {
  assert( regn_vlist.size()==0 && fR.size()==0 ) ;
  const int nRV = 7 ;
  int nR = 0;
  for ( int I=0; I<nx; ++I ) {
    for ( int J=0; J<ny; ++J ) {
      // set data of local L-region indexed by ir 
      for ( int ir=0; ir<2; ++ir ) {
	++nR ;
	regn_vlist.push_back( nRV ) ;
	for ( int il=0; il<nRV; ++il ) {
	  int iV = vrtx_index( ir, il, I, J ) ;
	  regn_vlist.push_back( iV ) ;
	}
	fR.push_back( UNSET ) ;
      }
    }
  }
  regn_vlist.push_back( nR ) ;
}

void build_lshape_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0, 
			double wx=0., double wy=0. ) {

  vector<int>    fV ;
  vector<double> xV, yV ;
  vector<int>    regn_vlist, fR ;
  
  LShapeMesh lshape_mesh(nx,ny,mesh_flag,wx,wy) ;
  lshape_mesh.build_vertex_coords  ( xV, yV, fV ) ;
  lshape_mesh.build_region_topology( regn_vlist, fR ) ;

  int nR = 2*nx*ny ;
  int nV = 5*nx*ny + 3*nx + ny +1 ;
  
  assert( xV.size()==nV ) ;
  assert( yV.size()==nV ) ;
  assert( fV.size()==nV ) ;
  assert( regn_vlist[regn_vlist.size()-1]==nR ) ;
  assert( fR.size()==nR ) ;
  
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR L-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("LShape-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
} ;

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// ROTATED MESH OF REGULAR TRIANGLES (standard flag: 5xx)
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
/*
     m2     m2+1    ...     m2+nd
      ________________________
     /\      /\      /\      /\
    /  \ D  /  \    /  \    /  \
   /    \  /    \  /    \  /    \
  /  C   \/      \/      \/      \
  --------------------------------
 m1     m1+1    ...    m1+nd    m1+nd+1
  \      /\      /\      /\      /
   \ B  /  \    /  \    /  \    /
    \  / A  \  /    \  /    \  /
     \/------\/------\/------\/
     m0      m0+1   ...      m0+nd  

     m1 = (m0+nd)+1
     m2 = (m1+nd+1)+1

     A = { m0,   m0+1, m1+1 } + i
     B = { m0,   m1+1, m1+0 } + i
     C = { m1,   m1+1, m2   } + i
     D = { m1+1, m2+1, m2   } + i

     (the vertices m1 and m1+nd+1 are taken as boundary nodes)
 */
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// input:
// ------
// nx, assigned to nd
// ny, wx, wy, mesh_flag, not used
//
void build_rotated_tria_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
			      double wx=0., double wy=0. ) {
  
  vector<double> xV, yV ;
  vector<int>    regn_vlist ;
  vector<int>    fV, fR ;

  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  int nd = nx ; // nd must be even
  double dx = 1./double(nd) ;
  double dy = dx ;
  for ( int j=0; j<=nd; ++j ) {
    for ( int i=0; i<=nd; ++i ) {
      xV.push_back( double(i)*dx ) ;
      yV.push_back( double(j)*dy ) ;
      fV.push_back( UNSET ) ;
    }
    if ( j<nd) {
      xV.push_back( 0. ) ;
      yV.push_back( double(j)*dy+dy/2. ) ;
      fV.push_back( UNSET ) ;
      for ( int i=1; i<=nd; ++i ) {
	xV.push_back( double(i)*dx-dx/2. ) ;
	yV.push_back( double(j)*dy+dy/2. ) ;
	fV.push_back( UNSET ) ;
      }
      xV.push_back( 1. ) ;
      yV.push_back( double(j)*dy+dy/2. ) ;
      fV.push_back( UNSET  ) ;
    }
  }

  int nR = 0 ;
  int m0 = 0 ;
  int m1 = (m0+nd)+1 ;
  int m2 = (m1+nd+1)+1 ;
  for ( int j=0; j<nd; ++j ) {
    for ( int i=0; i<nd; ++i ) {   // triangles of type A
      regn_vlist.push_back(3) ;
      regn_vlist.push_back(m0+0+i) ;
      regn_vlist.push_back(m0+1+i) ;
      regn_vlist.push_back(m1+1+i) ;
      fR.push_back(UNSET) ;
      ++nR ;
    }
    for ( int i=0; i<=nd; ++i ) {  // triangles of type B
      regn_vlist.push_back(3) ;
      regn_vlist.push_back(m0+0+i) ;
      regn_vlist.push_back(m1+1+i) ;
      regn_vlist.push_back(m1+0+i) ;
      fR.push_back(UNSET) ;
      ++nR ;
    }
    for ( int i=0; i<=nd; ++i ) {  // triangles of type C
      regn_vlist.push_back(3) ;
      regn_vlist.push_back(m1+0+i) ;
      regn_vlist.push_back(m1+1+i) ;
      regn_vlist.push_back(m2+0+i) ;
      fR.push_back(UNSET) ;
      ++nR ;
    }
    for ( int i=0; i<nd; ++i ) {  // triangles of type D
      regn_vlist.push_back(3) ;
      regn_vlist.push_back(m1+1+i) ;
      regn_vlist.push_back(m2+1+i) ;
      regn_vlist.push_back(m2+0+i) ;
      fR.push_back(UNSET) ;
      ++nR ;
    }
    m0 = m2 ;
    m1 = m0+(nd+1) ;
    m2 = m1+(nd+1)+1 ;
  }
  regn_vlist.push_back(nR) ;
  
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  int lstr(16) ;
  char str_nx[lstr] ; sprintf(str_nx,  "%d",nx) ;
  char str_ny[lstr] ; sprintf(str_ny,  "%d",ny) ;
  string mesh_name = string("Rotated-Mesh-500-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// DOUBLE QUADRILATERAL MESH (standard flag: 6xx)
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*

  nx=1 ==> nd=2

  i=0               i=2*np
  |----|----|----|----|
  |    |    |    |    |
  |----|----|----|----|
  |    |    |    |    |
  o----|----|----|----|
  |         |         |
  |         |         |
  |         |         |
  |---------|---------|
  i=0               i=np

  o --> nVB, (x0T,y0T)

  1) two different maps are used for bottom/top regions
  2) cmap_T (top) is taken from j=0 to j=nd, Y-coordinates are translated of 0.5

*/
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void build_double_quad_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
		             double wx=0., double wy=0. ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int np  = 2*nx ;                        // ny not used
  int nV  = (np+1)*(np+1)+3*nx*(2*nx+1) ; //
  int nR  = 10*nx*nx ;                    //
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  int flag = mesh_flag%10 ; 
  CoordinateMap cmap_B(  np,  np,flag,wx,wy) ; // bottom mesh
  CoordinateMap cmap_T(2*np,2*np,flag,wx,wy) ; // top    mesh
  int    nVB = cmap_B.vrtx_index(0,nx) ;
  double x0T = cmap_B.sxy(0,nx) ;
  double y0T = cmap_B.txy(0,nx) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;

  // bottom grid nodes
  for ( int j=0; j<=nx-1; ++j ) {
    for ( int i=0; i<=np; ++i ) {
      int iV = cmap_B.vrtx_index(i,j) ;
      xV[iV] = cmap_B.sxy(i,j) ;
      yV[iV] = cmap_B.txy(i,j) ;
      fV[iV] = UNSET ;
    }
  }
  { // interface nodes
    int j = 0 ;
    for ( int i=0; i<=2*np; i+=2 ) {
      int iV = cmap_T.vrtx_index(i,j) + nVB ;
      xV[iV] = cmap_B.sxy(i/2,nx) ;
      yV[iV] = 0.5 ; // (randomize) cmap_B.txy(i/2,nx) ;
      fV[iV] = UNSET ;
    }
    for ( int i=1; i<=2*np; i+=2 ) {
      int iV  = cmap_T.vrtx_index(i,  j) + nVB ;
      int iVm = cmap_T.vrtx_index(i-1,j) + nVB ;
      int iVp = cmap_T.vrtx_index(i+1,j) + nVB ;
      xV[iV] = 0.5*( xV[iVm]+xV[iVp] ) ;
      yV[iV] = 0.5*( yV[iVm]+yV[iVp] ) ;
      fV[iV] = UNSET ;
    }
  }
  // top grid nodes
  for ( int j=1; j<=2*nx; ++j ) {
    for ( int i=0; i<=2*np; ++i ) {
      int iV = cmap_T.vrtx_index(i,j) + nVB ;
      xV[iV] = j==2*nx ? double(i)/double(2*np) : cmap_T.sxy(i,j) + x0T ;
      yV[iV] = j==2*nx ?                     1. : cmap_T.txy(i,j) + y0T ;
      fV[iV] = UNSET ;
    }
  }
//-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR(nR) ;

  int iR = 0;
  for ( int j=0; j<nx-1; ++j ) {
    for ( int i=0; i<np; ++i ) {
      int iVa = cmap_B.vrtx_index(i,  j  ) ;
      int iVb = cmap_B.vrtx_index(i+1,j  ) ;
      int iVc = cmap_B.vrtx_index(i+1,j+1) ;
      int iVd = cmap_B.vrtx_index(i  ,j+1) ;
      regn_vlist.push_back(4) ;
      regn_vlist.push_back(iVa) ;
      regn_vlist.push_back(iVb) ;
      regn_vlist.push_back(iVc) ;
      regn_vlist.push_back(iVd) ;
      fR[iR++] = UNSET ;
    }
  }
  {
    int j = nx-1 ;
    for ( int i=0; i<np; ++i ) {
      int iVa = cmap_B.vrtx_index(i,  j  ) ;
      int iVb = cmap_B.vrtx_index(i+1,j  ) ;
      int iVc = cmap_T.vrtx_index(2*i+2,0) + nVB ;
      int iVd = cmap_T.vrtx_index(2*i+1,0) + nVB ;
      int iVe = cmap_T.vrtx_index(2*i  ,0) + nVB ;
      regn_vlist.push_back(5) ;
      regn_vlist.push_back(iVa) ;
      regn_vlist.push_back(iVb) ;
      regn_vlist.push_back(iVc) ;
      regn_vlist.push_back(iVd) ;
      regn_vlist.push_back(iVe) ;
      fR[iR++] = UNSET ;
    }
  }
  for ( int j=0; j<2*nx; ++j ) {
    for ( int i=0; i<2*np; ++i ) {
      int iVa = cmap_T.vrtx_index(i,  j  ) + nVB ;
      int iVb = cmap_T.vrtx_index(i+1,j  ) + nVB ;
      int iVc = cmap_T.vrtx_index(i+1,j+1) + nVB ;
      int iVd = cmap_T.vrtx_index(i  ,j+1) + nVB ;
      regn_vlist.push_back(4) ;
      regn_vlist.push_back(iVa) ;
      regn_vlist.push_back(iVb) ;
      regn_vlist.push_back(iVc) ;
      regn_vlist.push_back(iVd) ;
      fR[iR++] = UNSET ;
    }
  }
  regn_vlist.push_back(nR) ;

  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",600+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx)   ;
  string mesh_name = string("DoubleQuad-")+str_flag+string("-")+str_nx+string("x")+str_nx ;
  mesh . set_mesh_name( mesh_name ) ;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// Non-Convex Octagonal-BASED MESH (standard mflag: 8xy), x=[0|1], y=[0..3]
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class OctoMap {
private:
  int    mflag, nx, ny ;
  double dx, dy, wi, wj ;
  double sgn_i, sgn_j ;
public:
  OctoMap( int _nx=1, int _ny=1, int _mflag=0, double _wi=0., double _wj=0. ) : 
    nx(_nx), ny(_ny), mflag(_mflag), wi(_wi), wj(_wj) {
    sgn_i = +1. ;
    sgn_j = +1. ;
    dx    = 1./double(nx) ;
    dy    = 1./double(ny) ;
    srandom(10) ;
  }
  ~OctoMap() {}
  inline int  hor_face_index( int i, int j ) const { return i +  nx    * j ; }
  inline int  vrt_face_index( int i, int j ) const { return j +  ny    * i ; }
  inline int  vrtx_index    ( int i, int j ) const { return i + (nx+1) * j ; }
  inline bool is_bnd_vertex ( int i, int j ) const {
    return i==0 || i==nx || j==0 || j==ny ;
  }
  inline bool is_i_on_boundary( int i ) const { return i==0 || i==nx ; }
  inline bool is_j_on_boundary( int j ) const { return j==0 || j==ny ; }
  inline bool is_i_on_vrt_bndr( int i ) const { return i==nx/2 ;       }
  inline bool is_j_on_hor_bndr( int j ) const { return j==ny/2 ;       }
  inline double xy_rand( double a=-1., double b=1. ) const {
    //return a + (b-a) * double(random())/double(LONG_MAX) ;
    return a + (b-a) * double(random())/double(RAND_MAX) ;
  }
  // set sign for alternating shift
  void change_sign( int sflag ) {
    switch ( sflag )  {
    case 0: break ;
    case 1: sgn_i *= -1. ; break ;
    case 2: sgn_j *= -1. ; break ;
    default: assert(false) ;
    }
  }
  // map-0, identity
  inline double s0( double x, double y ) const { return x ; }
  inline double t0( double x, double y ) const { return y ; }
  // map-1, horizontal/vertical shift
  inline double s1( double x, double y ) const { return x+sgn_i*0.5*dx*wi ; }
  inline double t1( double x, double y ) const { return y+sgn_j*0.5*dy*wj ; }
  // map-2, (randomized vertices)
  inline double s2( double x, double y ) const { 
    return x + wi * xy_rand( -dx/2., dx/2. ) ;
  }
  inline double t2( double x, double y ) const { 
    return y + wj * xy_rand( -dy/2., dy/2. ) ;
  }
  // map-3, 
  inline double s3( double x, bool do_shift ) const { 
    return do_shift ? x+sgn_i*0.5*dx*wi : x ; 
  }
  inline double t3( double y, bool do_shift ) const { 
    return do_shift ? y+sgn_j*0.5*dy*wj : y ; 
  }

  inline double sxy( double x, double y, bool do_shift=false ) const {
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=s0(x,y) ; break ;
    case 1: retval=s1(x,y) ; break ;
    case 2: retval=s2(x,y) ; break ;
    case 3: retval=s3(x,do_shift) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
  inline double txy( double x, double y, bool do_shift=false ) const { 
    double retval(0.) ;
    switch( mflag ) {
    case 0: retval=t0(x,y) ; break ;
    case 1: retval=t1(x,y) ; break ;
    case 2: retval=t2(x,y) ; break ;
    case 3: retval=t3(y,do_shift) ; break ;
    default: assert(false) ;
    }
    return retval ;
  }
} ;

void build_NonConvexOcto_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, int mesh_flag=0,
			       double wx=0., double wy=0. ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  double  dx = 1./double(nx) ;
  double  dy = 1./double(ny) ;
  int  mflag = mesh_flag % 10 ;
  int  sflag = mesh_flag / 10 ;
  OctoMap       octo(nx,ny,mflag,wx/2.,wy/2.) ;
  CoordinateMap cmap(nx,ny,mflag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF OCTO-MESH
  //-----------------------------------------------------------------------------------------
  int nV = (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  // nodes of quads 
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      double xij = dx * double(i) ;
      double yij = dy * double(j) ;
      int iV = octo.vrtx_index(i,j) ;
      if ( octo.is_i_on_boundary(i) ||  octo.is_j_on_boundary(j) ){ 
	xV[iV] = xij ;
	yV[iV] = yij ; 
//    } else if ( octo.is_i_on_vrt_bndr(i) ) {
//      xV[iV] = xij ;
//      yV[iV] = cmap.txy( xij, yij ) ; // yij
      } else {
	xV[iV] = cmap.sxy( xij, yij ) ; // xij
	yV[iV] = cmap.txy( xij, yij ) ; // yij
      }
      fV[iV] = UNSET ;
    }
  }
  //LINE(--) ;
  int nV_hori = (nx+1)*(ny+1) ;
  // nodes of horizontal faces
  for ( int j=0; j<=ny; ++j ) {
    octo.change_sign( sflag ) ;
    for ( int i=0; i<nx; ++i ) {
      double xij = dx * ( double(i) + 0.5 ) ;
      double yij = dy *   double(j) ;
      int iV = nV_hori + octo.hor_face_index(i,j) ;
      if ( octo.is_j_on_boundary(j) ) {
	xV[iV] = xij ;
	yV[iV] = yij ; 
      } else {
	xV[iV] = octo.sxy( xij, yij ) ;
	yV[iV] = octo.txy( xij, yij, true ) ;
      }
      fV[iV] = UNSET ;
      //VAL(i) ; VAL(j) ; VAL(iV) ; VAL(xij) ; PRT(yij) ;
    }
  }
  //LINE(--) ;
  int nV_vert = nV_hori + nx*(ny+1) ;
  // nodes of vertical faces
  for ( int i=0; i<=nx; ++i ) {
    octo.change_sign( sflag ) ;
    for ( int j=0; j<ny; ++j ) {
      double xij = dx *   double(i) ;
      double yij = dy * ( double(j) + 0.5 ) ;
      int iV = nV_vert + octo.vrt_face_index(i,j) ;
      if ( octo.is_i_on_boundary(i) ) {
	xV[iV] = xij ;
	yV[iV] = yij ; 
//    } else if ( octo.is_i_on_vrt_bndr(i) ) {
//  	xV[iV] = xij ;
//  	yV[iV] = octo.txy( xij, yij ) ;
      } else {
	xV[iV] = octo.sxy( xij, yij, true ) ;
	yV[iV] = octo.txy( xij, yij ) ;
      }
      fV[iV] = UNSET ;
      //VAL(i) ; VAL(j) ; VAL(iV) ; VAL(xij) ; PRT(yij) ;
    }
  }
  //LINE(--) ;
  //for ( int i=0; i<fV.size(); ++i ) {
  //  VALA(i,xV) ; PRTA(i,yV) ; 
  //}
  //LINE(--) ;

  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF OCTO-MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {
      regn_vlist.push_back( 8 ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i,  j  )         ) ;
      regn_vlist.push_back( octo.hor_face_index(i,  j  )+nV_hori ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i+1,j  )         ) ;
      regn_vlist.push_back( octo.vrt_face_index(i+1,j  )+nV_vert ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i+1,j+1)         ) ;
      regn_vlist.push_back( octo.hor_face_index(i,  j+1)+nV_hori ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i,  j+1)         ) ;
      regn_vlist.push_back( octo.vrt_face_index(i,  j  )+nV_vert ) ;
      // ---
      fR.push_back( UNSET ) ;
    }
  }
  int nR = nx*ny ;
  regn_vlist.push_back( nR ) ;
  assert( nR == regn_vlist.size()/9 ) ;

  //PRT_ARR( regn_vlist ) ;

  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR OCTOGONAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",800+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("Non-Convex Octo-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}

void build_NonConvexNonRandomOcto_mesh( mesh_2Dv & mesh, int nx=1, int ny=1, 
                                        int mesh_flag=0, double wx=0., double wy=0. ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  double  dx = 1./double(nx) ;
  double  dy = 1./double(ny) ;
  int  mflag = mesh_flag % 10 ;
  int  sflag = mesh_flag / 10 ;
  OctoMap octo(nx,ny,mflag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF OCTO-MESH
  //-----------------------------------------------------------------------------------------
  int nV = (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  // nodes of quads 
  for ( int i=0; i<=nx; ++i ) {
    for ( int j=0; j<=ny; ++j ) {
      double xij = dx * double(i) ;
      double yij = dy * double(j) ;
      int iV = octo.vrtx_index(i,j) ;
      xV[iV] = xij ; // cmap.sxy( xij, yij ) ;
      yV[iV] = yij ; // cmap.txy( xij, yij ) ;
      fV[iV] = UNSET ;
    }
  }
  //LINE(--) ;
  int nV_hori = (nx+1)*(ny+1) ;
  // nodes of horizontal faces
  for ( int j=0; j<=ny; ++j ) {
    octo.change_sign( sflag ) ;
    for ( int i=0; i<nx; ++i ) {
      double xij = dx * ( double(i) + 0.5 ) ;
      double yij = dy *   double(j) ;
      int iV = nV_hori + octo.hor_face_index(i,j) ;
      if ( octo.is_j_on_boundary(j) ) {
	xV[iV] = xij ;
	yV[iV] = yij ; 
      } else {
	xV[iV] = octo.sxy( xij, yij ) ;
	yV[iV] = octo.txy( xij, yij, true ) ;
      }
      fV[iV] = UNSET ;
      //VAL(i) ; VAL(j) ; VAL(iV) ; VAL(xij) ; PRT(yij) ;
    }
  }
  //LINE(--) ;
  int nV_vert = nV_hori + nx*(ny+1) ;
  // nodes of vertical faces
  for ( int i=0; i<=nx; ++i ) {
    octo.change_sign( sflag ) ;
    for ( int j=0; j<ny; ++j ) {
      double xij = dx *   double(i) ;
      double yij = dy * ( double(j) + 0.5 ) ;
      int iV = nV_vert + octo.vrt_face_index(i,j) ;
      if ( octo.is_i_on_boundary(i) ) {
	xV[iV] = xij ;
	yV[iV] = yij ; 
      } else {
	xV[iV] = octo.sxy( xij, yij, true ) ;
	yV[iV] = octo.txy( xij, yij ) ;
      }
      fV[iV] = UNSET ;
      //VAL(i) ; VAL(j) ; VAL(iV) ; VAL(xij) ; PRT(yij) ;
    }
  }
  //LINE(--) ;
  //for ( int i=0; i<fV.size(); ++i ) {
  //  VALA(i,xV) ; PRTA(i,yV) ; 
  //}
  //LINE(--) ;

  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF OCTO-MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {
      regn_vlist.push_back( 8 ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i,  j  )         ) ;
      regn_vlist.push_back( octo.hor_face_index(i,  j  )+nV_hori ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i+1,j  )         ) ;
      regn_vlist.push_back( octo.vrt_face_index(i+1,j  )+nV_vert ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i+1,j+1)         ) ;
      regn_vlist.push_back( octo.hor_face_index(i,  j+1)+nV_hori ) ;
      // ---
      regn_vlist.push_back( octo.vrtx_index    (i,  j+1)         ) ;
      regn_vlist.push_back( octo.vrt_face_index(i,  j  )+nV_vert ) ;
      // ---
      fR.push_back( UNSET ) ;
    }
  }
  int nR = nx*ny ;
  regn_vlist.push_back( nR ) ;
  assert( nR == regn_vlist.size()/9 ) ;

  //PRT_ARR( regn_vlist ) ;

  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR OCTOGONAL-SHAPED MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  int lstr(16) ;
  char str_flag[lstr] ; sprintf(str_flag,"%d",800+mesh_flag) ;
  char str_nx  [lstr] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [lstr] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("Non-Convex Octo-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
}




#endif // end of _MESH_2D_BUILTIN_HH
