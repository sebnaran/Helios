// ===========================================================================================================================

extern "C" {
  void mof_( double * xyc,    // nodal coordinates of the QUAD cell (can be non-convex)
             int    & r_nmat, // number of materials (since this is a mixed cell, nmat>=2)
             double * cref,   // reference centroids
             double * f,      // volume fraction
             int    * ord,    // material ordering
             int    & oflag,  // flag for ordering of material ( if "0", automatic ordering is provided)
             double * grad,   // user specified initial guess of MoF interface normal
             int    & gflag,  // flag for initial guess for grad (if "0", guess is provided MoF method)
             double & eps1,   // tolerance for volume fraction cutting (e.g. ~1.d-10)
             double & eps2,   // tolerance for MoF angle optimization (e.g. ~1.d-5)
             double * cact,   // actual (reconstructed) centroids, e.g.
             double * cerr,   // centroid err, ||(x_ref - x_act)/len||^2
             double * xyt,    // coordinates of pure sub-triangles
             int    * nt,     // number of triangles for each materials
             double * xyp,    // nodal coordinates of polygonal sub-cell
             int    * np,     // number of polygon for each materials
             int    * nvp,    // number of vertices for each polygon
             int    & nit1,   // iteration counters
             int    & nit2,   // iteration counters
             int    & cyl,    // flag for selecting "x-y" [0] or "r-z" [2] coordinates
             double * xyg,    // global set of coordinates (where idxt() and idxp() is indexing)
             int    * idxp,   // index vector for pure polygons
             int    * idxt,   // index vector for pure triangles:
             double * ang,    // defines interface normal
             double * dst     // defines signed interface distance from origin
             ) ;
}

// "paint" a generic 2D curve on a mesh
class sub_mesh2Dv_painter_MOF {
private:
  static const int DIM  = 2 ;
  //static const int n_mat_intf = 2 ;
  static const int undef_zone = -1 ;
  const double epsi ;
  // const int nmat, nintf ;

  struct mof_node_struct {
  public:
    double xV, yV ;
    int    iV ;
    mof_node_struct( int _iV, double _xV, double _yV ) :
      iV(_iV), xV(_xV), yV(_yV) {}
    ~mof_node_struct() {}
    bool operator< ( const mof_node_struct & mof_node ) const { 
      bool retval = 
	xV < mof_node.xV || 
	( abs(xV-mof_node.xV) < 1e-14 && yV < mof_node.yV ) ;
      return retval ; 
    }
  } ;
  
  vector<SubGridCell> subgridcell_vec ;

  // output
  //  10 --> should be n_mat_intf
  // 128 --> max number of triangles
  FullMatrix grad, cact, xyt, xyp ;
  Vector cerr ;
  VecInt nt, np, nvp ;
  VecInt idxt, idxp, ord ;
  FullMatrix xyg ;
  Vector ang, dst ;

  // multiple lines
  vector<double> x0_vec, y0_vec, x1_vec, y1_vec ;
  
private:
  int    icrv ; // flag to select the curve
  double ax, ay, a0, tol ;
  double x0, y0, rd ;
  double x1 ;
  mesh_2Dv & mesh ;
  
  int nx, ny ;
  FullMatrix xv, yv ;
  Vector  vrtx_izone_vec ;
  VecBool mixed_cell ;
  vector<int> mixed_cell_vec ;
  
  // ---------------------------------------
  // ---------------------------------------
  
  // single curve
  int  vrtx_izone( int iV ) ;
  int  vrtx_izone( double xV, double yV ) ;
  void build_vrtx_izone_vec() ;
  void set_subcell_structure() ;
  // ---
  double curve( double x, double y ) ;
  double line     ( double x, double y ) ;
  double circle   ( double x, double y ) ;
  double ellipse  ( double x, double y ) ;
  double hyperbola( double x, double y ) ;
  
  // strip
  int  strip_izone( int iV ) ;
  int  strip_izone( double x, double y ) ;
  void build_vrtx_izone_strip_vec() ;
  void set_subcell_structure_strip() ;
  // ---
  double vertical_strip( int flag, double x, double y ) ;
  double circular_strip( int flag, double x, double y ) ;

  // -----------
  int  get_izone_flag( double x, double y ) ;
  void get_mixed_cell_list() ;
  void process_mixed_cell_list() ;
  void compute_volume_fraction( int iR, vector<int> & regn_vlist ) ;
  void compute_local_subgrid( int iR, vector<int> & regn_vlist ) ;
  void split_the_cell( int iR, Vector & f, FullMatrix & cref, FullMatrix & xyc ) ;
  void filter_mof_nodes( int nnode, FullMatrix & xyg, SubGridCell & subc, VecInt & vrtx_index_map ) ;
  void filter_mof_nodes( SubGridCell & subc ) ;
  void set_node_flag   ( SubGridCell & subc ) ;
  
  void setup() ;

public:
  sub_mesh2Dv_painter_MOF( mesh_2Dv & _mesh, int _n_mat_intf=1 ) : mesh(_mesh), epsi(1e-12), nx(100), ny(100), 
								   icrv(-1), n_mat_intf(_n_mat_intf)
  {
    setup() ;
  }
  ~sub_mesh2Dv_painter_MOF() {}

  // single curve
  void paint_a_line      ( double _ax, double _ay, double _a0, double _tol ) ;
  void paint_a_circle    ( double _x0, double _y0, double _rd              ) ;
  void paint_an_ellipse  ( double _x0, double _y0, double _ax, double _ay  ) ;
  void paint_an_hyperbola( double _x0, double _y0, double _ax, double _ay  ) ;

  // strips
  void paint_vertical_strip( double x0=0.45, double y0=0.55 ) ;
  void paint_circular_strip() ;
  
  // multiple lines
  void paint_multiple_lines( vector<double> & _x0_vec, vector<double> & _y0_vec, 
			     vector<double> & _x1_vec, vector<double> & _y1_vec ) ;
  void set_subcell_structure_mlines() ;
  void build_vrtx_izone_mlines_vec() ;
  int  mlines_izone( int iV ) ;
  int  mlines_izone( double xV, double yV ) ;
  double interface_line( int imat, double x, double y ) ;

  // number of subcells
  int n_subcell() { return subgridcell_vec.size() ; }
  SubGridCell & get_subcell( int ic ) {
    assert( 0<=ic && ic<subgridcell_vec.size() ) ;
    return subgridcell_vec[ic] ; 
  }
} ;

void sub_mesh2Dv_painter_MOF::setup() {
  xv.setup(nx+1,ny+1) ;
  yv.setup(nx+1,ny+1) ;
  // ---
  int nR = mesh.n_region() ;
  mixed_cell.setup( nR ) ;
  // ---
  mixed_cell_vec.resize(0) ;
  // ---
  int nV = mesh.n_vertex() ;
  vrtx_izone_vec.setup( nV ) ;

  // output
  //  10 --> should be n_mat_intf
  // 128 --> max number of triangles
  grad.setup(DIM,nmat) ;
  cact.setup(DIM,nmat) ;
  xyt .setup(DIM,128) ;
  xyp .setup(DIM,128) ;
  cerr.setup(nmat) ;
  nt  .setup(10) ;
  np  .setup(10) ;
  nvp .setup(10*2) ;
  idxt.setup(128) ;
  idxp.setup(128) ; 
  ord .setup(10)  ;
  xyg .setup(DIM,128) ;
  ang .setup(nmat) ;
  dst .setup(nmat) ;

  // check quads
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = mesh.n_regn_face(iR) ;
    assert( nRF==4 ) ;
  }
}

void sub_mesh2Dv_painter_MOF::get_mixed_cell_list() {
  //MSGF("begin  sub_mesh2Dv_painter_MOF::get_mixed_cell_list") ;

  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    
    // check if MOF can be used
    int nRF = mesh.n_regn_face(iR) ;
    assert( nRF==4 ) ;

    // get vertex list
    vector<int> regn_vlist ;
    mesh.get_regn_vrtx( iR, regn_vlist ) ;
    int nRV = regn_vlist.size() ;
    
    // set mixed_flag
    // get the first vrtx_flag different than undef_zone
    mixed_cell(iR) = false ;
    int vrtx_flag = undef_zone ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV  = regn_vlist[ilV] ; 
      int izn = vrtx_izone_vec( regn_vlist[ilV] ) ;
      if ( izn!=undef_zone ) {
	vrtx_flag = izn ;
	break ;
      }
    }

    // check if iR is a mixed region
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = regn_vlist[ilV] ; 
      if ( vrtx_izone_vec(iV)!=vrtx_flag && vrtx_izone_vec(iV)!=undef_zone ) {
	mixed_cell_vec.push_back(iR) ;
	mixed_cell(iR) = true ;
	break ;
      }
    }
  }
#if 0 // debug TBR
  PRT(icrv) ;
  PRT_ARR(mixed_cell_vec) ;
  exit(0) ;
#endif
  //MSGF("end of sub_mesh2Dv_painter_MOF::get_mixed_cell_list") ;
}

// compute MOF's input and call it
// this routine works only for two materials
void sub_mesh2Dv_painter_MOF::process_mixed_cell_list() {
  //MSGF("begin  sub_mesh2Dv_painter_MOF::process_mixed_cell_list") ;

  int n_mixed_cells = mixed_cell_vec.size() ;
  for ( int il=0; il<n_mixed_cells; ++il ) {

    // get region's id
    int iR = mixed_cell_vec[il] ;
    
    // get vertex list
    vector<int> regn_vlist ;
    mesh.get_regn_vrtx( iR, regn_vlist ) ;

    // build local subgrid 
    compute_local_subgrid( iR, regn_vlist ) ; 

    // compute volume fractions and center coordinates
    PRT(nmat) ;
    vector<double> area_vec(nmat), xc_vec(nmat), yc_vec(nmat) ;
    //TBR
    //double area_0(0.), area_1(0.) ; 
    //double xc_0(0.), yc_0(0.), xc_1(0.), yc_1(0.) ;
    for ( int s=0; s<nx; ++s ) {
      for ( int t=0; t<ny; ++t ) {
	double xc  = ( xv(s,t) + xv(s+1,t) + xv(s+1,t+1) + xv(s,t+1) )/4. ;
	double yc  = ( yv(s,t) + yv(s+1,t) + yv(s+1,t+1) + yv(s,t+1) )/4. ;
	double dx0 = xv(s+1,t) - xv(s,t) ;
	double dy0 = yv(s+1,t) - yv(s,t) ;
	// ---
	double dx1 = xv(s,t+1) - xv(s,t) ;
	double dy1 = yv(s,t+1) - yv(s,t) ;
	// ---
	double area = abs( dx0*dy1 - dx1*dy0 ) ;
	// ---
	int izn_c = get_izone_flag(xc,yc)>0 ; 
	// ---
	area_vec[izn_c] += area ;
	xc_vec  [izn_c] += area * xc ;
	yc_vec  [izn_c] += area * yc ;
	
	// TBR
	// 	if ( izn_c==0 ) {
	// 	  area_0 += area ; 
	// 	  xc_0   += area * xc ;
	// 	  yc_0   += area * yc ;
	// 	} else if ( izn_c==1 ) { 
	// 	  area_1 += area ; 
	// 	  xc_1   += area * xc ;
	// 	  yc_1   += area * yc ;
	// 	} else {
	// 	  // not implemented yet
	// 	  assert(false) ;
	// 	}
      }
    }
    // ----
    int nv = 4 ;
    FullMatrix xyc(DIM,nv) ;
    // ---
    xyc(0,0) = xv(0,0) ;
    xyc(1,0) = yv(0,0) ;
    // ---
    xyc(0,1) = xv(nx,0) ;
    xyc(1,1) = yv(nx,0) ;
    // ---
    xyc(0,2) = xv(nx,ny) ;
    xyc(1,2) = yv(nx,ny) ;
    // ---
    xyc(0,3) = xv(0,ny) ;
    xyc(1,3) = yv(0,ny) ;
    // ----
    FullMatrix cref(DIM,nmat) ;
    // ----
    for ( int imat=0; imat<nmat; ++imat ) {
      cref(0,imat) = xc_vec[imat] / area_vec[imat] ;
      cref(1,imat) = yc_vec[imat] / area_vec[imat] ;
    }
    //     cref(0,0) = xc_0 / area_0 ; // TBR
    //     cref(1,0 )= yc_0 / area_0 ;
    //     // ----
    //     cref(0,1) = xc_1 / area_1 ;
    //     cref(1,1) = yc_1 / area_1 ;
    //     // ----
    double mR = mesh.get_regn_measure( iR ) ;
    //double diff = 1. - ( area_0 + area_1 )/mR ; // TBR
    Vector f(nmat) ;
    double sum_area(0.) ;
    for ( int imat=0; imat<nmat-1; ++imat ) {
      f(imat) = area_vec[imat]/mR ;
      sum_area += area_vec[imat]/mR ;
    }
    f(nmat-1) = 1. - sum_area ;
    //f(0) = area_0/mR + diff/2. ;
    //f(1) = 1. - f(0) ;
    printf("iR=%3i f0=%14.7e  f1=%14.7e\n",iR,f(0),f(1)) ;
    // --- call MOF
    double tol = 1e-2 ;
    //double tol = 1e-10 ;
    if ( f(0)>=tol && f(1)>=tol ) {
      split_the_cell( iR, f, cref, xyc ) ;
    }
  }

  //MSGF("end of sub_mesh2Dv_painter_MOF::process_mixed_cell_list") ;
}

// remove replicated nodes
void sub_mesh2Dv_painter_MOF::filter_mof_nodes( int nnode, FullMatrix & xyg, SubGridCell & subc, VecInt & vrtx_index_map ) {
  //MSGF("begin sub_mesh2Dv_painter_MOF::filter_mof_nodes nmat = "<<nmat) ;
  
  double dist_tol = 1e-6 ;

  //MSGF("begin  filter_mof_nodes") ;
  
  vector<mof_node_struct> node_list ;
  for ( int iV=0; iV<nnode; ++iV ) {
    node_list.push_back( mof_node_struct( iV, xyg(0,iV),  xyg(1,iV) ) ) ;
  }

  sort( node_list.begin(), node_list.end() ) ;

  int ip = 0 ;
  vrtx_index_map.setup(nnode) ;

  { // first item
    int i=0 ;
    subc.xV.push_back( node_list[i].xV ) ;
    subc.yV.push_back( node_list[i].yV ) ;
    vrtx_index_map(node_list[i].iV) = ip ;    
  }
  for ( int i=1; i<nnode; ++i ) {
    double dx = node_list[i].xV - subc.xV[ip] ;
    double dy = node_list[i].yV - subc.yV[ip] ;
    double dist = sqrt( pow(dx,2) + pow(dy,2) ) ;
    if ( dist > dist_tol ) {
      subc.xV.push_back( node_list[i].xV ) ;
      subc.yV.push_back( node_list[i].yV ) ;

      ++ip ;
      vrtx_index_map(node_list[i].iV) = ip ;    
    } else {
      vrtx_index_map(node_list[i].iV) = ip ;    
    }
  }

  //MSGF("end of sub_mesh2Dv_painter_MOF::filter_mof_nodes") ;
}

// remove the hanging nodes 
// (hanging nodes are present from output of MOF)
// note that fR is never modified
void sub_mesh2Dv_painter_MOF::filter_mof_nodes( SubGridCell & subc ) {
  //MSGF("begin sub_mesh2Dv_painter_MOF::filter_mof_nodes") ;

  int nV = subc.xV.size() ;
  VecInt vrtx_remap(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { vrtx_remap(iV) = -1 ; }

  vector<double> new_xV, new_yV ;
  vector<int>    new_fV ;
  int new_nV = 0 ;
  int ip = 0 ;
  int nR = subc.regn_vlist.back() ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRV = subc.regn_vlist[ip++] ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = subc.regn_vlist[ip++] ;
      if ( vrtx_remap(iV) == -1 ) { 
	vrtx_remap(iV) = new_nV++ ;
	new_xV.push_back( subc.xV[iV] ) ;
	new_yV.push_back( subc.yV[iV] ) ;
	new_fV.push_back( subc.fV[iV] ) ;
      }
    }
  }
  assert( new_nV<=nV ) ;
  
  // debug, TBR
  //for ( int iV=0; iV<nV; ++iV ) {
  //  printf("iV=%2i new_iV=%2i\n",iV,vrtx_remap(iV)) ;
  //}
 
  vector<int> new_regn_vlist ;
  ip = 0 ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRV = subc.regn_vlist[ip++] ;
    new_regn_vlist.push_back( nRV ) ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = subc.regn_vlist[ip++] ;
      new_regn_vlist.push_back( vrtx_remap(iV) ) ;
    }
  }
  new_regn_vlist.push_back( nR ) ;

  assert( new_regn_vlist.size()==subc.regn_vlist.size() ) ;

  // debug, TBR
  //PRT_ARR(subc.regn_vlist) ;
  //LINE(--) ;
  //PRT_ARR(new_regn_vlist) ;
  //exit(0) ;

  subc.regn_vlist.resize(0) ;
  subc.xV.resize(0) ;
  subc.yV.resize(0) ;
  subc.fV.resize(0) ;
  
  for ( int i=0; i<new_regn_vlist.size(); ++i ) {
    subc.regn_vlist.push_back( new_regn_vlist[i] ) ;
  }
  
  for ( int i=0; i<new_xV.size(); ++i ) {
    subc.xV.push_back( new_xV[i] ) ;
    subc.yV.push_back( new_yV[i] ) ;
    subc.fV.push_back( new_fV[i] ) ;
  }

  //MSGF("end of sub_mesh2Dv_painter_MOF::filter_mof_nodes") ;
}

void sub_mesh2Dv_painter_MOF::set_node_flag( SubGridCell & subc ) {
  int iR = subc.cell_id ;

  // get vertex list
  vector<int> regn_vlist ;
  mesh.get_regn_vrtx( iR, regn_vlist ) ;
  int nRV = regn_vlist.size() ;

  // get face list
  vector<int> regn_flist ;
  mesh.get_regn_face( iR, regn_flist ) ;
  int nRF = regn_flist.size() ;

  // ...
  for ( int il=0; il<subc.xV.size(); ++il ) {
    double xl = subc.xV[il] ;
    double yl = subc.yV[il] ;
    bool ok_vrtx = false ;
    // check if il is a base_mesh vertex
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int    iV = regn_vlist[ilV] ;
      double xV = mesh.coords_V(iV,0) ;
      double yV = mesh.coords_V(iV,1) ;
      if ( !ok_vrtx && abs(xl-xV) + abs(yl-yV)<1e-12 ) {
	subc.fV[il] = -iV-1 ;
	ok_vrtx = true ;
      }
    }
    // if not, check the base_mesh face where it is located
    if ( !ok_vrtx ) {
      
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	int    iV0 = regn_vlist[ilV] ;
	int    iV1 = regn_vlist[ (ilV+1)%nRV ] ;
	// ---
	double x0 = mesh.coords_V(iV0,0) ;
	double y0 = mesh.coords_V(iV0,1) ;
	double x1 = mesh.coords_V(iV1,0) ;
	double y1 = mesh.coords_V(iV1,1) ;
	// ---
	double ddV = (y0-yl)*(xl-x1) - (x0-xl)*(yl-y1) ;
	if ( abs(ddV)<1.e-12 ) {
	  subc.fV[il] = regn_flist[ilV] ;
	}
      }
    }
  }

#if 0 // debug TBR
  LINE(--) ;
  PRT(iR) ;
  for ( int il=0; il<subc.xV.size(); ++il ) {
    PRTA(il,subc.fV);
  }
  exit(0) ;
#endif
}

void sub_mesh2Dv_painter_MOF::split_the_cell( int iR, Vector & f, FullMatrix & cref, FullMatrix & xyc ) {  
  //MSGF("begin  sub_mesh2Dv_painter_MOF::split_the_cell -->> iR = "<<iR) ; 

  // input data
  double eps1 = 1e-10 ;
  double eps2 = 1e-5  ;
  
  int gflag = 0 ;
  int oflag = 0 ;

  int nit1  = 0 ;
  int nit2  = 0 ;
  int cyl   = 0 ;
  
  // calling MOF subroutine
  //MSGF("calling MOF subroutine...") ;
  int r_nmat = nmat ;
  mof_( xyc.add(),    // nodal coordinates of the QUAD cell (can be non-convex)
        r_nmat,       // number of materials (since this is a mixed cell, nmat>=2)
        cref.add(),   // reference centroids
        f.add(),      // volume fraction
        ord.add(),    // material ordering
        oflag,        // flag for ordering of material ( if "0", automatic ordering is provided)
        grad.add(),   // user specified initial guess of MoF interface normal
        gflag,        // flag for initial guess for grad (if "0", guess is provided MoF method)
        eps1,         // tolerance for volume fraction cutting (e.g. ~1.d-10)
        eps2,          // tolerance for MoF angle optimization (e.g. ~1.d-5)
        cact.add(),   // actual (reconstructed) centroids, e.g.
        cerr.add(),   // centroid err, ||(x_ref - x_act)/len||^2
        xyt.add(),    // coordinates of pure sub-triangles
        nt.add(),     // number of triangles for each materials
        xyp.add(),    // nodal coordinates of polygonal sub-cell
        np.add(),     // number of polygon for each materials
        nvp.add(),    // number of vertices for each polygon
        nit1,         // iteration counters
        nit2,         // iteration counters
        cyl,          // flag for selecting "x-y" [0] or "r-z" [2] coordinates
        xyg.add(),    // global set of coordinates (where idxt() and idxp() is indexing)
        idxp.add(),   // index vector for pure polygons
        idxt.add(),   // index vector for pure triangles:
        ang.add(),    // defines interface normal
        dst.add()     // defines signed interface distance from origin
        ) ;
  //MSGF("out from MOF subroutine") ;

  int ntria = 0 ;
  for ( int i=0; i<nmat; ++i ) {
    ntria += nt(i) ;
    //PRTV(i,nt) ;
  }
  int nnode = 3*ntria ;

  int npoly      = 0 ;
  int nnode_poly = 0 ;
  for ( int i=0; i<nmat; ++i ) {
    npoly += np(i) ;
    for ( int ip=0; ip<np(i); ++ip ) {
      nnode_poly += nvp(ip) ;
      //PRTV(ip,nvp) ;
    }
  }

  subgridcell_vec.push_back( SubGridCell(iR) ) ;
  SubGridCell & subc = subgridcell_vec.back() ;

  VecInt vrtx_index_map ;

  // remove duplicated nodes
  filter_mof_nodes( nnode, xyg, subc, vrtx_index_map ) ;

  for ( int i=0; i<subc.xV.size(); ++i ) { subc.fV.push_back(0) ; }

  int offset = 0 ;
  int nptot  = 0 ;
  int nR     = 0 ;
  for ( int k=0; k<nmat; ++k ) {

    if ( np(k)==1 ) {

      nR++ ;
      subc.regn_vlist.push_back( nvp(nptot) ) ;
      subc.fR.push_back( k ) ;
      for ( int i=0; i<nvp(nptot); ++i ) {
        int iV = vrtx_index_map( idxp(offset)-1 ) ;
        subc.regn_vlist.push_back( iV ) ;
        offset++ ;
      }

    } else {

      nR++ ;
      subc.regn_vlist.push_back( nvp(nptot) ) ;
      subc.fR.push_back( k ) ;
      for ( int i=0; i<nvp(nptot); ++i ) {
        int iV = vrtx_index_map( idxp(offset)-1 ) ;
        subc.regn_vlist.push_back( iV ) ;
        offset++ ;
      }

      nR++ ;
      subc.regn_vlist.push_back( nvp(nptot) ) ;
      subc.fR.push_back( k ) ;
      for ( int i=0; i<nvp(nptot+1); ++i ) {
        int iV = vrtx_index_map( idxp(offset)-1 ) ;
        subc.regn_vlist.push_back( iV ) ;
        offset++ ;
      }

    }

    nptot += np(k) ;

  }
  subc.regn_vlist.push_back( nR ) ;

  // remove hanging nodes
  filter_mof_nodes( subc ) ;

  // reset fV
  set_node_flag( subc ) ;

  //assert( subc.regn_vlist.back()==2 ) ;
  //assert( subc.regn_vlist[0]==4 ) ;
  //assert( subc.regn_vlist[5]==4 ) ;

  if (0) {
    int nV = subc.xV.size() ;
    PRT(nV) ;
    for ( int iV=0; iV<nV; ++iV ) {
      printf("iV=%2i xV=%14.7e yV=%14.7e fV=%2i\n",iV,subc.xV[iV],subc.yV[iV],subc.fV[iV]) ;
    }
    int nR = subc.regn_vlist.back() ;
    PRT(nR) ;
    PRT_ARR(subc.regn_vlist) ;
    exit(0) ;
  }

  //MSGF("end of sub_mesh2Dv_painter_MOF::split_the_cell") ; 
}

// compute temporary subgrid
void sub_mesh2Dv_painter_MOF::compute_local_subgrid( int iR, vector<int> & regn_vlist ) {  
  // get local coordinates of cell iR
  int iVA = regn_vlist[0] ;
  int iVB = regn_vlist[1] ;
  int iVC = regn_vlist[2] ;
  int iVD = regn_vlist[3] ;

  double xA = mesh.coords_V( iVA, 0 ) ; 
  double xB = mesh.coords_V( iVB, 0 ) ; 
  double xC = mesh.coords_V( iVC, 0 ) ; 
  double xD = mesh.coords_V( iVD, 0 ) ; 

  double yA = mesh.coords_V( iVA, 1 ) ; 
  double yB = mesh.coords_V( iVB, 1 ) ; 
  double yC = mesh.coords_V( iVC, 1 ) ; 
  double yD = mesh.coords_V( iVD, 1 ) ;

  // build grid
  // ---
  double ds = 1./double(nx) ;
  double dt = 1./double(ny) ;
  for ( int s=0; s<=nx; ++s ) {
    double xs_AB = xA + ds*double(s) * ( xB-xA ) ;
    double ys_AB = yA + ds*double(s) * ( yB-yA ) ;
    // ----
    double xs_DC = xD + ds*double(s) * ( xC-xD ) ;
    double ys_DC = yD + ds*double(s) * ( yC-yD ) ;
    // ----
    for ( int t=0; t<=ny; ++t ) {
      xv(s,t) = xs_AB + dt*double(t)*( xs_DC-xs_AB ) ;
      yv(s,t) = ys_AB + dt*double(t)*( ys_DC-ys_AB ) ;
    }
  }
}

// =================================================================
// =================================================================

int sub_mesh2Dv_painter_MOF::get_izone_flag( double x, double y ) {
  int retval = -1 ;
  if      (  0<=icrv && icrv<=3  ) { return vrtx_izone  (x,y) ; }
  else if (  4<=icrv && icrv<=5  ) { return strip_izone (x,y) ; }
  else if ( 10<=icrv && icrv<=20 ) { return mlines_izone(x,y) ; }
  return retval ;
}

// =================================================================
// =================================================================
//
// single curve
//
// =================================================================
// =================================================================

int sub_mesh2Dv_painter_MOF::vrtx_izone( int iV ) {
  assert( 0<=iV && iV<mesh.n_vertex() ) ;
  double xV = mesh.coords_V( iV, 0 ) ;
  double yV = mesh.coords_V( iV, 1 ) ;
  return vrtx_izone( xV, yV ) ;
}

int sub_mesh2Dv_painter_MOF::vrtx_izone( double xV, double yV ) {
  double rxy = curve( xV, yV ) ;
  int retval = undef_zone ;         // undef  zone, probable intersection
  if ( rxy<-epsi ) { retval = 0 ; } // first  zone
  if ( rxy>+epsi ) { retval = 1 ; } // second zone
  return retval ;
}

void sub_mesh2Dv_painter_MOF::build_vrtx_izone_vec() {
  int nV = mesh.n_vertex() ;
  assert( vrtx_izone_vec.size()==nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    vrtx_izone_vec(iV) = vrtx_izone(iV) ;
  }
}

void sub_mesh2Dv_painter_MOF::set_subcell_structure() {
  // set the zone vertex flag
  build_vrtx_izone_vec() ;

  // find which cells are mixed
  get_mixed_cell_list() ;

  // process the mixed cells
  process_mixed_cell_list() ;
}

double sub_mesh2Dv_painter_MOF::curve( double x, double y ) {
  double retval(0.) ;
  switch( icrv ) {
  case 0: retval = line     (x,y) ; break ; 
  case 1: retval = circle   (x,y) ; break ; 
  case 2: retval = ellipse  (x,y) ; break ; 
  case 3: retval = hyperbola(x,y) ; break ; 
  default: assert(false) ;
  }
  return retval ;
}

double sub_mesh2Dv_painter_MOF::line( double x, double y ) {
  return ax*x + ay*y + a0 ;
}

double sub_mesh2Dv_painter_MOF::circle( double x, double y ) {
  return pow(x0-x,2) + pow(y0-y,2) - pow(rd,2 ) ;
}

double sub_mesh2Dv_painter_MOF::ellipse( double x, double y ) {
  return pow( (x0-x)/ax, 2 ) + pow( (y0-y)/ay, 2 ) - 1. ;
}

double sub_mesh2Dv_painter_MOF::hyperbola( double x, double y ) {
  return pow( (x0-x)/ax, 2 ) - pow( (y0-y)/ay, 2 ) - 1. ;
}

void sub_mesh2Dv_painter_MOF::paint_a_line( double _ax, double _ay, double _a0, double _tol ) {
  MSGF("begin  sub_mesh2Dv_painter_MOF::paint_a_line") ;
  // set line coefficients
  ax  = _ax ;
  ay  = _ay ;
  a0  = _a0 ;
  tol = _tol ;
  // ---
  assert( ax>0. || abs(ax)<epsi && ay>0. ) ;
  // ---
  icrv = 0 ;
  // ---
  set_subcell_structure() ;
  MSGF("end of sub_mesh2Dv_painter_MOF::paint_a_line") ;
}

void sub_mesh2Dv_painter_MOF::paint_a_circle( double _x0, double _y0, double _rd ) {
  // set circle's center and radius
  x0 = _x0 ;
  y0 = _y0 ; 
  rd = _rd ;
  // ---
  icrv = 1 ;
  // ---
  set_subcell_structure() ;
}

void sub_mesh2Dv_painter_MOF::paint_an_ellipse( double _x0, double _y0, double _ax, double _ay ) {
  // set ellipse's center and axis
  x0 = _x0 ;
  y0 = _y0 ; 
  ax = _ax ;
  ay = _ay ;
  // ---
  icrv = 2 ;
  // ---
  set_subcell_structure() ;
}

void sub_mesh2Dv_painter_MOF::paint_an_hyperbola( double _x0, double _y0, double _ax, double _ay ) {
  // set hyperbola coeffcients
  x0 = _x0 ;
  y0 = _y0 ; 
  ax = _ax ;
  ay = _ay ;
  // ---
  icrv = 3 ;
  // ---
  set_subcell_structure() ;
}

// =================================================================
// =================================================================
//
// strip
//
// =================================================================
// =================================================================

// include set_subcell_structure
void sub_mesh2Dv_painter_MOF::paint_vertical_strip( double _x0, double _x1 ) {
  MSGF("begin  sub_mesh2Dv_painter_MOF::paint_vertical_strip") ;
  icrv = 4 ;
  x0 = _x0 ;
  x1 = _x1 ;
  set_subcell_structure_strip() ;
  MSGF("end of sub_mesh2Dv_painter_MOF::paint_vertical_strip") ;
}

void sub_mesh2Dv_painter_MOF::paint_circular_strip() {
  MSGF("begin  sub_mesh2Dv_painter_MOF::paint_circular_strip") ;
  icrv = 5 ;
  set_subcell_structure_strip() ;
  MSGF("end of sub_mesh2Dv_painter_MOF::paint_circular_strip") ;
}

void sub_mesh2Dv_painter_MOF::set_subcell_structure_strip() {
  // set the zone vertex flag
  build_vrtx_izone_strip_vec() ;
    
  // find which cells are mixed
  get_mixed_cell_list() ; 

  // process the mixed cells
  process_mixed_cell_list() ;
}

void sub_mesh2Dv_painter_MOF::build_vrtx_izone_strip_vec() {
  int nV = mesh.n_vertex() ;
  assert( vrtx_izone_vec.size()==nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    vrtx_izone_vec(iV) = strip_izone(iV) ;
   }
}

int sub_mesh2Dv_painter_MOF::strip_izone( int iV ) {
  assert( 0<=iV && iV<mesh.n_vertex() ) ;
  double xV = mesh.coords_V( iV, 0 ) ;
  double yV = mesh.coords_V( iV, 1 ) ;
  return strip_izone( xV, yV ) ;
}

int sub_mesh2Dv_painter_MOF::strip_izone( double xV, double yV ) {
  double r0(0.), r1(0.) ;
  if ( icrv==4 ) {
    r0 = vertical_strip( 0, xV, yV ) ;
    r1 = vertical_strip( 1, xV, yV ) ;
  } else if ( icrv==5 ) {
    r0 = circular_strip( 0, xV, yV ) ;
    r1 = circular_strip( 1, xV, yV ) ;
  }
  // ------------------------------------
  int retval = undef_zone ;                  // undef  zone, probable intersection
  if ( 0.<=r0 && r1<=epsi ) { retval = 1 ; } // second zone
  else                      { retval = 0 ; } // first  zone
  // ------------------------------------
  return retval ;
}

// (same code in Test_Dirichlet_1)
double sub_mesh2Dv_painter_MOF::vertical_strip( int flag, double x, double y ) {
  //double x0(0.45), x1(0.55) ;
  double retval(0.) ;
  if      ( flag==0 ) { retval = x-x0 ; }
  else if ( flag==1 ) { retval = x-x1 ; }
  return retval ;
}

// (same code in Test_Dirichlet_1)
double sub_mesh2Dv_painter_MOF::circular_strip( int flag, double x, double y ) {
  double xc(-0.5), yc(0.5)  ;
  double x0(0.45), x1(0.55) ;
  // ------------------------
  double retval(0.) ;
  if ( flag==0 ) {
    double r0 = sqrt( pow(x0-xc,2) + pow(yc,2) ) ;
    retval = pow(xc-x,2) + pow(yc-y,2) - pow(r0,2) ;
  } else if ( flag==1 ) {
    double r1 = sqrt( pow(x1-xc,2) + pow(yc,2) ) ;
    retval = pow(xc-x,2) + pow(yc-y,2) - pow(r1,2) ;
  }
  // ------------------------
  return retval ;
}

// =================================================================


// =================================================================
// =================================================================
//
// multiple lines
//
// =================================================================
// =================================================================

// include set_subcell_structure
void sub_mesh2Dv_painter_MOF::paint_multiple_lines( vector<double> & _x0_vec, vector<double> & _y0_vec, 
						    vector<double> & _x1_vec, vector<double> & _y1_vec ) {
  MSGF("begin  sub_mesh2Dv_painter_MOF::paint_multiple_lines") ;

  PRT(n_mat_intf) ;
  PRT(_x0_vec.size()) ;

  assert( n_mat_intf == _x0_vec.size() ) ;
  for ( int imat=0; imat<n_mat_intf; ++imat ) {
    x0_vec.push_back( _x0_vec[imat] ) ;
    y0_vec.push_back( _y0_vec[imat] ) ;
    x1_vec.push_back( _x1_vec[imat] ) ;
    y1_vec.push_back( _y1_vec[imat] ) ;
  }

  icrv = 10 ;
  set_subcell_structure_mlines() ;
  MSGF("end of sub_mesh2Dv_painter_MOF::paint_multiple_lines") ;
}

// always using mlines
void sub_mesh2Dv_painter_MOF::set_subcell_structure_mlines() {
  MSGF("begin  sub_mesh2Dv_painter_MOF::set_subcell_structure_mlines") ;
  // set the zone vertex flag
  build_vrtx_izone_mlines_vec() ;

  // find which cells are mixed
  get_mixed_cell_list() ; 

  // process the mixed cells
  process_mixed_cell_list() ;
  MSGF("end of sub_mesh2Dv_painter_MOF::set_subcell_structure_mlines") ;
}

void sub_mesh2Dv_painter_MOF::build_vrtx_izone_mlines_vec() {
  //MSGF("begin  sub_mesh2Dv_painter_MOF::build_vrtx_izone_mlines_vec") ;
  int nV = mesh.n_vertex() ;
  assert( vrtx_izone_vec.size()==nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    vrtx_izone_vec(iV) = mlines_izone(iV) ;
    VAL(iV) ; PRTV(iV,vrtx_izone_vec) ;
  }
  //MSGF("end of sub_mesh2Dv_painter_MOF::build_vrtx_izone_mlines_vec") ;
}

int sub_mesh2Dv_painter_MOF::mlines_izone( int iV ) {
  //MSGF("begin  sub_mesh2Dv_painter_MOF::mlines_izone -->> "<<iV) ;
  assert( 0<=iV && iV<mesh.n_vertex() ) ;
  double xV = mesh.coords_V( iV, 0 ) ;
  double yV = mesh.coords_V( iV, 1 ) ;
  //MSGF("end of sub_mesh2Dv_painter_MOF::mlines_izone") ;
  return mlines_izone( xV, yV ) ;
}

int sub_mesh2Dv_painter_MOF::mlines_izone( double xV, double yV ) {
  //MSGF("begin  sub_mesh2Dv_painter_MOF::mlines_izone -->> "<<xV<<"  "<<yV) ;
  assert( icrv==10 ) ;
  int retval = n_mat_intf ; // undef  zone, probable intersection
  for ( int imat=0; imat<n_mat_intf; ++imat ) {
    double ri = interface_line( imat, xV, yV ) ;
    if ( ri<0. ) { retval = imat ; break ; }
  }
  //MSGF("end of sub_mesh2Dv_painter_MOF::mlines_izone") ;
  return retval ;
}

double sub_mesh2Dv_painter_MOF::interface_line( int imat, double x, double y ) {
  assert( 0<=imat && imat<x0_vec.size() ) ;
  double x0 = x0_vec[imat] ;
  double y0 = y0_vec[imat] ;
  double x1 = x1_vec[imat] ;
  double y1 = y1_vec[imat] ;
  double rv = ( x-x0 )*( y1-y ) - ( x1-x )*( y-y0 ) ;

  //MSG("interface_line -->> ") ; VAL(x) ; VAL(y) ; PRT(rv) ;

  return rv ;
}

// =================================================================
