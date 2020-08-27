#ifndef _MESH2DC_HH
#define _MESH2DC_HH

// PATCH FOR 2D REGION WITH CURVED FACES
const int  MAX_N_NODE    = 4 ;
const bool CURVED_FACES  = true ;
const bool NO_VERTICAL   = false ;
const bool NO_HORIZONTAL = false ;

// Pi = X(iV0) + alpha_i |F| tng_F + delta_i |F| nor_F
// X(t) = sum_{i=0}^{N} Pi phi_i(t),  t \in [0,1]
// phi_i( alpha_j ) = \delta_{ij}, i,j=0...N
class curved_face {
private:
  vector<double> alpha_vec ;
  vector<double> delta_vec ;
public:
  curved_face() {}
  ~curved_face() {}
  void set_n_node(int new_n_node)  { 
    assert( 0<new_n_node && new_n_node<=MAX_N_NODE ) ;
    alpha_vec.resize(new_n_node) ;
    delta_vec.resize(new_n_node) ;
  }
  void check_input ( int i ) const {
    assert( 0<=i && i<n_node() ) ;
  }
  int n_node() const { 
    assert( alpha_vec.size()==delta_vec.size() ) ; 
    return alpha_vec.size() ; 
  }
  double   alpha(int i) const { check_input(i) ; return alpha_vec[i] ; }
  double   delta(int i) const { check_input(i) ; return delta_vec[i] ; }
  double & alpha(int i)       { check_input(i) ; return alpha_vec[i] ; }
  double & delta(int i)       { check_input(i) ; return delta_vec[i] ; }
  bool is_curved_face() { return alpha_vec.size()>2 ; } 
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
  dataMat cvd_R_coords ;
  dataMat cvd_R_area   ;
  // ---

private:
  vector<curved_face> cvd_flist ;

private:
  // return a random real between a and b (a<b)
  double xy_random( double a=-1., double b=1. ) {
    assert( a<b ) ;
    return a + (b-a) * double(random())/double(INT_MAX) ;
  }

  // aux method
  double phi( int iN, int nFN ) { 
    assert( 0<=iN && iN<n_curved_faces() ) ;
    return iN==0 || iN==nFN-1 ? 0. : 1. ; 
  }

  // initialization methods
  void setup_cvd_flist() ;
  void setup_region_geometry() ;
  void setup_region_geometry( int iR ) ;

  // setup geometry of curved faces
  void check_region_geometry() ;
  void check_region_geometry( int iR ) ;
  
  // local nodes giving the shape of a polynomial face
  double face_node( int iN, int iF, int idim ) ;

private:
  // thes quantities are for a global check of consistency
  double max_diff_area, max_diff_bary_X, max_diff_bary_Y ;

public:
  mesh_2Dc() {}
  ~mesh_2Dc() {}

  // used in main to init curved faces
  void init_curve() ;

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
  
  // area of region R  (possibly with curved faces)
  double get_cvd_regn_measure( int iR ) ;
} ;

// used in main to init curved faces
void mesh_2Dc :: init_curve() {
  setup_cvd_flist() ;
  setup_region_geometry() ;
  check_region_geometry() ;
}

// every face has two extrema, so TWO points is the minimum
bool mesh_2Dc :: is_curved_face( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  return cvd_flist[iF].n_node() > 2 ;
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

// position vector of iN-th node of face iF
// node iN is given by the pair ( alpha(iN), delta(iN) )
double mesh_2Dc :: face_node( int iN, int iF, int idim ) {
  assert( 0<=iF && iF<cvd_flist.size() ) ;
  assert( 0<=iN && iN<cvd_flist[iF].n_node() ) ;
  //assert( is_curved_face(iF) ) ;
  int    iV0  = face_vrtx( iF, 0 ) ;
  double lenF = get_face_measure( iF ) ;
  return coords_V( iV0, idim ) 
    + lenF * get_tng( iF, idim ) * cvd_flist[iF].alpha(iN) 
    + lenF * get_nor( iF, idim ) * cvd_flist[iF].delta(iN) ;
}

// parameterization of the curved face F through X(t)
// t in [0,1] --> X(t) on F
// X(t)     = sum_{0<=i<N} X_i phi_i(t)
// X_i      = face_node( iN, iF, idim) 
// phi_i(t) = prod_{j!=i} (t-alpha_j)/(alpha_i-alpha_j)
// phi_i( alpha_j ) = 1 if i==j; =0 otherwise
// all loops are from iN=0 (first) to iN=nFN-1 (last)
double mesh_2Dc :: face_pnts( double t, int iF, int idim ) {
  //assert( 0<=t    && t<=1        ) ;
  assert( 0<=iF   && iF<n_face() ) ;
  assert( 0<=idim && idim<DIM    ) ;

  //if (!is_curved_face(iF) ) {
  //  PRT(iF) ;
  //  assert( is_curved_face(iF)     ) ;
  //}
  double retval =0. ;
  int nFN = cvd_flist[iF].n_node() ;
  for ( int iN=0; iN<nFN; ++iN ) {                 // summation loop on retval
    double phi_iN = 1. ;
    for ( int jN=0; jN<nFN; ++jN ) {               // product loop on phi_iN
      double alpha_i = cvd_flist[iF].alpha(iN) ;
      double alpha_j = cvd_flist[iF].alpha(jN) ;
      if ( iN!=jN ) {
	phi_iN *= (t-alpha_j)/(alpha_i-alpha_j) ;
	// 	phi_iN *= 
	// 	  ( t                      -cvd_flist[iF].alpha(jN) ) / 
	// 	  ( cvd_flist[iF].alpha(iN)-cvd_flist[iF].alpha(jN) ) ;
      }
    }
    // accumulate: retval <-- X_iN * phi_iN(t)
    retval += face_node( iN, iF, idim ) * phi_iN ;
  }
  return retval ;
}

// distance of the projection of X(t) onto the face from the first point
double mesh_2Dc :: face_alpha( double t, int iF ) {
  assert( 0<=t  && t<=1        ) ;
  assert( 0<=iF && iF<n_face() ) ;
  return t ;
}

// distance of X(t) from the planar face through iN==0 (first) and iN==nFN-1 (last)
double mesh_2Dc :: face_delta( double t, int iF ) {
  //assert( 0<=t    && t<=1        ) ;
  assert( 0<=iF   && iF<n_face() ) ;
  assert( is_curved_face(iF)     ) ;
  double retval =0. ;
  int nFN = cvd_flist[iF].n_node() ;
  for ( int iN=1; iN<nFN-1; ++iN ) {               // summation loop on retval
    double phi_iN = 1. ;
    for ( int jN=0; jN<nFN; ++jN ) {               // product loop on phi_iN
      double alpha_i = cvd_flist[iF].alpha(iN) ;
      double alpha_j = cvd_flist[iF].alpha(jN) ;
      if ( iN!=jN ) {
	phi_iN *= (t-alpha_j)/(alpha_i-alpha_j) ;
      }
    }
    // accumulate: retval <-- delta * phi_iN(t)
    retval += cvd_flist[iF].delta(iN) * phi_iN ;
  }
  return retval ;
}

// X'(t) \in F, t \in [0,1]
// phi_i (t) = PROD_{j!=i} (t-alpha_j)/(alpha_i-alpha_j)
// phi_i'(t) = num_i(t) / den_i 
// num_i = SUM_{k!=i} PROD_{j!=i,k} (t-alpha_j)
// den_i = PROD_{j!=i} 1/(alpha_i-alpha_j)
double mesh_2Dc :: face_derv( double t, int iF, int idim ) {
  assert( 0<=iF   && iF<n_face() ) ;
  assert( 0<=idim && idim<DIM    ) ;
  //assert( is_curved_face(iF)     ) ;
  
  double retval = 0. ;
  int nFN = cvd_flist[iF].n_node() ;
  for ( int iN=0; iN<nFN; ++iN ) {
    
    // compute denominator of dphi_i
    double dphi_den = 1. ;
    for ( int jN=0; jN<nFN; ++jN ) {
      if ( iN!=jN ) {
	dphi_den *= cvd_flist[iF].alpha(iN)-cvd_flist[iF].alpha(jN) ;
      }
    }
    
    // compute numerator of dphi_i
    double dphi_num = 0. ;
    for ( int kN=0; kN<nFN; ++kN ) {
      if ( kN!=iN ) {
	double prod =1. ;
	for ( int jN=0; jN<nFN; ++jN ) {
	  if ( jN!=iN && jN!=kN ) {
	    prod *= t-cvd_flist[iF].alpha(jN) ;
	  }
	}
	dphi_num += prod ;
      }
    }
    
    // accumulate X'(t) = sum_i X_i dphi_i(t)
    retval += face_node( iN, iF, idim ) * dphi_num / dphi_den ;
  }
  //if ( iF==3 && idim==1 ) exit(0) ;
  return retval ;
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

// through this setup every internal face may be curved
// for every mesh face there is one entry in cvd_flist
void mesh_2Dc :: setup_cvd_flist() {
  cvd_flist.resize(nF) ;
  for ( int iF=0; iF<nF; ++iF ) {

    // set face length
    double lenF = get_face_measure(iF) ;

    // set number of face's node
    int nFN = 2 ;    
    if ( is_internal_face(iF) && CURVED_FACES ) { // CURVED FACES
      // very simple check to avoid some pathologies of a triangular mesh
      bool has_tria_neigh = bool(false) ;
      bool has_bnd_vrtx   = bool(false) ;
      for ( int s=0; s<2; ++s ) {
	int iR = face_regn( iF, s ) ;
	has_tria_neigh |= n_regn_face(iR)==3 ;
	int iV = face_vrtx( iF, s ) ;
	has_bnd_vrtx |= is_boundary_vrtx( iV ) ;
      }
      
      // only vertical faces may be curve
      bool is_vert_face = abs( coords_V( face_vrtx(iF,0), 0 ) -
			       coords_V( face_vrtx(iF,1), 0 ) )<1.e-14 ;

      bool is_horz_face = abs( coords_V( face_vrtx(iF,0), 1 ) -
			       coords_V( face_vrtx(iF,1), 1 ) )<1.e-14 ;

      bool check_curved = ( !has_tria_neigh || !has_bnd_vrtx ) ;
      check_curved &= NO_VERTICAL   ? !is_vert_face : true ;
      check_curved &= NO_HORIZONTAL ? !is_horz_face : true ;
      
      // only a face that does not have a triangular neighbour 
      // and is not on the boundary can be set as curved using nFN nodes
      if ( check_curved ) {
	nFN = min( 2 + int(random() % (MAX_N_NODE-1)), MAX_N_NODE ) ;
	//nFN = MAX_N_NODE ;
      }
    } // end of --> if ( is_internal_face(iF) ) {...

    // set data of face parameterization
    curved_face & F = cvd_flist[iF] ;
    F.set_n_node(nFN) ;
    // --------------
    int den = nFN/2 ;
    for ( int iN=0; iN<nFN; ++iN ) {
      F.alpha(iN) = double(iN)/double(nFN-1) ;
      if ( nFN<=3 ) { F.delta(iN) = 0.125 * double( iN*(nFN-1-iN)/den ) ; }              // symmetric
      else          { F.delta(iN) = 0.125 * double( iN*(nFN-1-iN)/den ) * double(iN) ; } // asymmetric
    }
  }
}

void mesh_2Dc :: setup_region_geometry() {
  MSG("begin setup_curved_regn_geom"<<endl<<flush) ;
  cvd_R_area  .setup(nR,1  ) ;
  cvd_R_coords.setup(nR,DIM) ;
  for ( int iR=0; iR<nR; ++iR ) {
    setup_region_geometry( iR ) ;
  }
  MSG("end-->setup_curved_regn_geom"<<endl<<flush) ;
}

void mesh_2Dc :: setup_region_geometry( int iR ) {
  assert( 0<=iR && iR<n_region() ) ;

  vector<int> regn_flist ;
  get_regn_face( iR, regn_flist ) ;
  const int nRF = n_regn_face( iR ) ;
  const int nRV = nRF ;
  const int nFV = DIM ;

  Gauss4 quad ;
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
  assert( abs(area_X-area_Y)<1.e-14 ) ;
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
  LINE(---) ;
  VAL( area_tot ) ; PRT( abs(1.-area_tot) ) ;
  VAL( xbary    ) ; PRT( abs(0.5-xbary)    ) ;
  VAL( ybary    ) ; PRT( abs(0.5-ybary)    ) ;
  LINE(---) ;
  PRT(max_diff_area) ; 
  PRT(max_diff_bary_X) ; 
  PRT(max_diff_bary_Y) ; 
  LINE(---) ;
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

      int n_face_node =  cvd_flist[iF].n_node() ;
      int max_np      = 2 * n_face_node ; // 2 is a safe factor
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

#endif // end of _MESH2DC_HH
