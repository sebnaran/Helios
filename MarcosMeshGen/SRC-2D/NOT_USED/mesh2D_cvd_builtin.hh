#ifndef _MESH_2D_CVD_BUILTIN_HH
#define _MESH_2D_CVD_BUILTIN_HH

#define CURVED_FACES 0

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// QUADRILATERAL-BASED MESH with CURVED faces (standard flag: 7xx)
//-----------------------------------------------------------------------------------------
// MESH FLAGS
// 700 --> default, straight edges
// 701 --> Los Alamos map, curvilinear edges
// 702 --> My map, curvilinear edges
// ---
// (standard coefficients follows)
//-----------------------------------------------------------------------------------------
void build_quad_cvd_mesh( mesh_2Dc & mesh, int nx=1, int ny=1, int mesh_flag=0,
			  double wx=0., double wy=0., int ilev=0 ) {
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  // cmap_flag takes the "x" in 70x
  int cmap_flag = mesh_flag % 10 ;
  CurvedMap cmap(nx,ny,cmap_flag,wx,wy) ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD CVD MESH
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
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (underlying skeleton)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  char str_flag[10] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nx  [10] ; sprintf(str_nx,  "%d",nx)   ;
  char str_ny  [10] ; sprintf(str_ny,  "%d",ny)   ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nx+string("x")+str_ny ;
  mesh . set_mesh_name( mesh_name ) ;
  //-----------------------------------------------------------------------------------------
  // set up curved face stuff
  //-----------------------------------------------------------------------------------------
#if CURVED_FACES
  // we assume that the boundary faces of a square are always straight segments
  mesh.setup_CurvedMap ( nx, ny, cmap_flag, wx, wy ) ;
  mesh.setup_cvd_flist() ;

  int    gl(-1), iF(-1) ;
  double s0(0.), s1(0.), gv(0.) ;
  double dx = cmap.get_dx() ;
  double dy = cmap.get_dy() ;

  const int HORIZONTAL_LINE = 0 ;
  const int VERTICAL_LINE   = 1 ;

  int iR = 0 ;
  for ( int j=0; j<ny; ++j ) {
    for ( int i=0; i<nx; ++i ) {

      vector<int> regn_flist ;
      mesh.get_regn_face( iR, regn_flist ) ;
      int nRF = regn_flist.size() ;
      assert( nRF==4 ) ;

      // face 0, (i,  j  )-->(i+1,j  ), s0=dx*i, s1=dx*(i+1), HORIZ
      iF = regn_flist[0] ;
      if ( mesh.face_regn( iF, 0 ) == iR && mesh.is_internal_face(iF) ) {
	s0 = dx*double(i)   ;
	s1 = dx*double(i+1) ;
	gl = HORIZONTAL_LINE ; // j
	gv = dy*double(j) ;
	mesh.setup_cvd_flist( iF, gl, gv, s0, s1 ) ;
      }

      // face 1, (i+1,j  )-->(i+1,j+1), s0=dy*j, s1=dy*(j+1), VERT
      iF = regn_flist[1] ;
      if ( mesh.face_regn( iF, 0 ) == iR && mesh.is_internal_face(iF) ) {
	s0 = dy*double(j)   ;
	s1 = dy*double(j+1) ;
	gl = VERTICAL_LINE ; // i+1 + n_HorizFaces ;
	gv = dx*double(i+1) ; 
	mesh.setup_cvd_flist( iF, gl, gv, s0, s1 ) ;
      }
      
      // face 2, (i+1,j+1)-->(i,j+1), s0=dx*i, s1=dx*(i+1), HORIZ
      iF = regn_flist[2] ;
      if ( mesh.face_regn( iF, 0 ) == iR && mesh.is_internal_face(iF) ) {
	s0 = dx*double(i+1)   ;
	s1 = dx*double(i) ;
	gl = HORIZONTAL_LINE ; // j+1 ;
	gv = dy*double(j+1) ;
	mesh.setup_cvd_flist( iF, gl, gv, s0, s1 ) ;
      }

      // face 3, (i,  j+1)-->(i+1,j), s0=dy*j, s1=dy*(j+1), VERT
      iF = regn_flist[3] ;
      if ( mesh.face_regn( iF, 0 ) == iR && mesh.is_internal_face(iF) ) {
	s0 = dy*double(j+1) ;
	s1 = dy*double(j) ;
	gl = VERTICAL_LINE ; // i + n_HorizFaces ;
	gv = dy*double(i) ;
	mesh.setup_cvd_flist( iF, gl, gv, s0, s1 ) ;
      }

      // new region
      ++iR ;
    }
  }
  //mesh.check_mesh() ;
#endif
}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// MULTILAYERED DISK WITH CURVED AZIMUTAL FACES
//-----------------------------------------------------------------------------------------
// MESH FLAGS used to select coaxial structures
// 710 --> default, single layer
// 711 --> two layers structure
// 712 --> real coaxial cable
// ---
// nx, wx, wy not used for the moment
// ny is used as multiplying factor in the internal refinement, set by mesh2D_init_cvd.hh
//-----------------------------------------------------------------------------------------
void set_vertex_ring( int i, CoaxMap & cmap,
		      vector<double> & xV, vector<double> & yV, vector<int> & fV ) {
  assert( xV.size()==yV.size() ) ;
  assert( xV.size()==fV.size() ) ;
  int nV = cmap.get_nV( i ) ;
  for ( int j=0; j<nV; ++j ) {
    int iV = cmap.vrtx_index(i,j) ;
    xV[iV] = cmap.xij(i,j) ;
    yV[iV] = cmap.yij(i,j) ;
    fV[iV] = UNSET ;
//     cout << " iV = " << setw(3)  << iV 
//     	 << " i  = " << setw(3)  << i  
//     	 << " j  = " << setw(3)  << j 
//     	 << " xV = " << setw(14) << setprecision(8) << xV[iV] 
//     	 << " yV = " << setw(14) << setprecision(8) << yV[iV] 
//     	 << endl ;
  }
}

void set_tria_cells( int igl, CoaxMap & cmap, 
		    vector<int> & regn_vlist, vector<int> & fR, int & nR ) { 
  assert( igl==0 ) ;
  //MSG("grid line --> ") ; PRT(igl) ;
  int ncell = cmap.get_nV(1) ;
  int nV1   = cmap.get_nV(1) ;
  int nRV = 3 ;
  for ( int j=0; j<ncell; ++j ) {
    ++nR ; // update cell's counter
    regn_vlist.push_back( nRV ) ;
    regn_vlist.push_back( cmap.vrtx_index( 0, 0 ) ) ;
    regn_vlist.push_back( cmap.vrtx_index( 1, j ) ) ;
    regn_vlist.push_back( cmap.vrtx_index( 1, (j+1)%nV1 ) ) ;
    fR.push_back( UNSET ) ;
    //MSG("  cell --> ") ; PRT(j) ;
    //cout << "   (0,0) --> (1," << j << ") --> (1," << (j+1)%nV1 << ") --> (0,0) " << endl ;
  }
}

void set_quad_cells( int igl, CoaxMap & cmap, 
		     vector<int> & regn_vlist, vector<int> & fR, int & nR ) { 
  //MSG("grid line --> ") ; PRT(igl) ;

  int nV0   = cmap.get_nV(igl)   ; // number of nodes on grid line igl
  int nV1   = cmap.get_nV(igl+1) ; // number of nodes on grid line igl+1
  int ncell = min( nV0, nV1 ) ;    // number of cells
  int nF0   = nV0 / ncell ;        // number of low   faces
  int nF1   = nV1 / ncell ;        // number of upper faces

  //VAL(nV0) ; PRT(nV1) ;
  //VAL(nF0) ; PRT(nF1) ;

  int i0 = igl   ;  
  int i1 = igl+1 ;  
  for ( int jc=0; jc<ncell; ++jc ) {
    int j0 = jc * nF0 ;
    int j1 = jc * nF1 ;

    // ---
    //MSG("  cell --> ") ; PRT(jc) ;
    // ---
    
    // update cell's counter
    ++nR ;

    // set number of cell nodes
    int nRV = nF0 + nF1 + 2 ; // low side (nF0+1) + upper side (nF1+1)
    regn_vlist.push_back( nRV ) ;

    // starting node, always (i0,j0)
    regn_vlist.push_back( cmap.vrtx_index( i0, j0 ) ) ;
    //cout << "   (" << i0 << "," << j0 << ")" ;
    // grid line i+1, all but the upper-last corner node
    for ( int k=0; k<nF1; ++k ) {
      regn_vlist.push_back( cmap.vrtx_index( i1, (j1+k)%nV1 ) ) ;
      //cout << " --> (" << i1 << "," << (j1+k)%nV1 << ")" ;
    }

    // the upper-last corner node 
    regn_vlist.push_back( cmap.vrtx_index( i1, (j1+nF1)%nV1 ) ) ;
    //cout << " --> (" << i1 << "," << (j1+nF1)%nV1 << ")" ;

    // grid line i, all but the starting node, only if igl>0
    for ( int k=0; k<nF0; ++k ) {
      regn_vlist.push_back( cmap.vrtx_index( i0, (j0+nF0-k)%nV0 ) ) ;
      //cout << " --> (" << i0 << "," << (j0+nF0-k)%nV0 << ")" ;
    }
    //cout << " --> (" << i0 << "," << j0 << ")" << endl ;
    fR.push_back( UNSET ) ;
  }
}

void set_cvd_tria_cells( int igl, CoaxMap & cmap, mesh_2Dc & mesh, int & iR ) { 
  const int RADIAL_LINE  = CoaxMap::RADIAL_LINE  ;
  const int AZIMUTH_LINE = CoaxMap::AZIMUTH_LINE ;

  assert( igl==0 ) ;
  //MSG("grid line, tria --> ") ; PRT(igl) ;

  int iF(-1), gl(-1) ;
  double r0(-1.), r1(-1.), t0(-1.), t1(-1.), gv(-1.) ;

  int ncell = cmap.get_nV(1) ;
  int nV1   = cmap.get_nV(1) ;
  int nRV = 3 ;
  for ( int j=0; j<ncell; ++j ) {
    
    vector<int> regn_flist ;
    mesh.get_regn_face( iR, regn_flist ) ;
    assert( regn_flist.size()==3 ) ;

#if 0 // RADIAL FACES ARE ALWAYS STRAIGHT LINES
    // (0,0)  -->(1,j),   radial
    iF = regn_flist[0] ;
    if ( mesh.face_regn( iF, 0 ) == iR ) {
      r0 = cmap.rho(0,0) ;
      r1 = cmap.rho(1,j) ;
      gl = RADIAL_LINE ; // j ;
      gv = double(j) * cmap.get_dt(1) ;
      mesh.setup_cvd_flist( iF, gl, gv, r0, r1 ) ;
      //VAL(iF) ; VAL(gl) ; VAL(gv) ; VAL(r0) ; PRT(r1) ;
    }
#endif

    // (1,j)  -->(1,j+1), circular
    iF = regn_flist[1] ;
    if ( mesh.face_regn( iF, 0 ) == iR ) {
      t0 = cmap.tht(1,j  ) ;
      t1 = cmap.tht(1,j+1) ;
      gl = AZIMUTH_LINE ;
      gv = cmap.rho(1,j) ;
      mesh.setup_cvd_flist( iF, gl, gv, t0, t1 ) ;
      //VAL(iF) ; VAL(gl) ; VAL(gv) ; VAL(t0) ; PRT(t1) ;
    }
    
#if 0 // RADIAL FACES ARE ALWAYS STRAIGHT LINES
    // (1,j+1)-->(0,0),   radial
    iF = regn_flist[2] ;
    if ( mesh.face_regn( iF, 0 ) == iR ) {
      r0 = cmap.rho(1,j+1) ;
      r1 = cmap.rho(0,0) ;
      gl = RADIAL_LINE ; // (j+1)%nV1 ;
      gv = double( (j+1)%nV1 ) * cmap.get_dt(1) ;
      mesh.setup_cvd_flist( iF, gl, gv, r0, r1 ) ;
      //VAL(iF) ; VAL(gl) ; VAL(gv) ; VAL(r0) ; PRT(r1) ;
    }
#endif
    
    // 
    ++iR ;
  }
}

void set_cvd_quad_cells( int igl, CoaxMap & cmap, mesh_2Dc & mesh, int & iR, 
			 bool ok_igl0=true, bool ok_igl1=true ) {
  const int RADIAL_LINE  = CoaxMap::RADIAL_LINE  ;
  const int AZIMUTH_LINE = CoaxMap::AZIMUTH_LINE ;

  assert( igl>0 ) ;
  //LINE(--) ;
  //MSG("grid line, quad --> ") ; PRT(igl) ;

  int iF(-1), gl(-1) ;
  double r0(-1.), r1(-1.), t0(-1.), t1(-1.), gv(-1.) ;

  int nV0   = cmap.get_nV(igl)   ; // number of nodes on grid line igl
  int nV1   = cmap.get_nV(igl+1) ; // number of nodes on grid line igl+1
  int ncell = min( nV0, nV1 ) ;    // number of cells
  int nF0   = nV0 / ncell ;        // number of low   faces
  int nF1   = nV1 / ncell ;        // number of upper faces

  //PRT(ncell) ;
  //VAL(nV0) ; PRT(nV1) ;
  //VAL(nF0) ; PRT(nF1) ;
  
  // jc : local cell index
  // iR : global cell index
  int i0 = igl   ;  
  int i1 = igl+1 ;  
  for ( int jc=0; jc<ncell; ++jc ) {
    int j0 = jc * nF0 ;
    int j1 = jc * nF1 ;

    // set the list of the region's faces
    vector<int> regn_flist ;
    mesh.get_regn_face( iR, regn_flist ) ;
    
    // set local face counter
    int ilF = 0 ;

    // (i,j)  -->  (i+1,j), radial
    iF = regn_flist[ ilF++ ] ;
    assert( mesh.face_regn(iF,0)==iR || mesh.face_regn(iF,1)==iR ) ;
#if 0 // RADIAL FACES ARE ALWAYS STRAIGTH LINES
    if ( mesh.face_regn( iF, 0 ) == iR ) {
      r0 = cmap.rho(i0,j0) ;
      r1 = cmap.rho(i1,j0) ;
      gl = RADIAL_LINE ; // j ;
      gv = double(j0) * cmap.get_dt( igl ) ;
      mesh.setup_cvd_flist( iF, gl, gv, r0, r1 ) ;
      //VAL(iF) ; VAL(iR) ; VAL(gl) ; VAL(gv) ; VAL(r0) ; PRT(r1) ;
    }
#endif
    
    // (i+1,j)  -->  (i+1,j+1), azimuth
    for ( int k=0; k<nF1; ++k ) {
      iF = regn_flist[ ilF++ ] ;
      assert( mesh.face_regn(iF,0)==iR || mesh.face_regn(iF,1)==iR ) ;
      if ( mesh.face_regn( iF, 0 ) == iR && ok_igl1 ) {
	t0 = cmap.tht( i1, (j1+k)   ) ;
	t1 = cmap.tht( i1, (j1+k+1) ) ;
	gl = AZIMUTH_LINE ;
	gv = cmap.rho( i1, j1 ) ;
	mesh.setup_cvd_flist( iF, gl, gv, t0, t1 ) ;
	//VAL(iF) ; PRT(iR) ; MSG("==>") ; VAL(k) ; VAL(gl) ; VAL(gv) ; VAL(t0) ; PRT(t1) ;
      }
    }

    // (i+1,j+1)  -->  (i,j+1), radial
    iF = regn_flist[ ilF++ ] ;
    assert( mesh.face_regn(iF,0)==iR || mesh.face_regn(iF,1)==iR ) ;
#if 0  // RADIAL FACES ARE ALWAYS STRAIGTH LINES
    if ( mesh.face_regn( iF, 0 ) == iR ) {
      r0 = cmap.rho(i1,j1) ;
      r1 = cmap.rho(i0,j1) ;
      gl = RADIAL_LINE ; // (j+1)%nt ;
      gv = double( j1+nF1 ) * cmap.get_dt( igl+1 ) ;
      mesh.setup_cvd_flist( iF, gl, gv, r0, r1 ) ;
      //VAL(iF) ; VAL(iR) ; VAL(gl) ; VAL(gv) ; VAL(r0) ; PRT(r1) ;
    }
#endif

    // (i,j+1)  -->  (i,j), azimuth
    for ( int k=0; k<nF0; ++k ) {
      iF = regn_flist[ ilF++ ] ;
      assert( mesh.face_regn(iF,0)==iR || mesh.face_regn(iF,1)==iR ) ;
      if ( mesh.face_regn( iF, 0 ) == iR && ok_igl0 ) {
	t0 = cmap.tht( i0, (j0+nF0-k) ) ;
	t1 = cmap.tht( i0, (j0+nF0-k-1) ) ;
	gl = AZIMUTH_LINE ;
	gv = cmap.rho( i0, j0 ) ;
	mesh.setup_cvd_flist( iF, gl, gv, t0, t1 ) ;
	//VAL(iF) ; PRT(iR) ; MSG("-->") ; VAL(k) ; VAL(gl) ; VAL(gv) ; VAL(t0) ; PRT(t1) ;
      }
    }

    // update cell counter
    ++iR ;
  }
}

// WARNING: nlev instead of ny, input as multiplying factor for refinements in coax 
//                              structure, set by mesh2D_init_cvd.hh
void build_disk_cvd_mesh( mesh_2Dc & mesh, int nx=1, int ny=1, int mesh_flag=0,
			  double wx=0., double wy=0., int ilev=0 ) {
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  // SET COAXIAL CABLE STRUCTURE
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  // cmap_flag takes the "x" in 71x
  int cmap_flag = mesh_flag % 10 ;
  int nlev = int( pow(2.,ilev) ) ;
  CoaxStructure coax( cmap_flag, nlev ) ;
  //PRT( coax.n_layer() ) ;
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  // SET COORDINATE MAP
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  CoaxMap cmap( coax ) ;
  //cmap.print_coax_grid() ; // debugging: print coaxial grid 
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF QUAD CVD MESH
  //-----------------------------------------------------------------------------------------
  int nV = cmap.n_vertex() ; 
  //PRT(nV) ;
  vector<int>    fV(nV) ;
  vector<double> xV(nV), yV(nV) ;
  // set origin

  int igl = 0 ;  
  int nl = cmap.n_layer() ;

  { // first layer
    int il = 0 ; // il=0 !!!

    set_vertex_ring( igl, cmap, xV, yV, fV ) ;
    ++igl ;
    
    // set other grid lines
    int nk = il==nl-1 ? cmap.get_nr(il) + 1 : cmap.get_nr(il) ;
    for ( int k=1; k<nk; ++k ) {
      set_vertex_ring( igl, cmap, xV, yV, fV ) ;
      ++igl ;
    }
  }

  // set other vertex layers
  for ( int il=1; il<nl; ++il ) {
    int nk = il==nl-1 ? cmap.get_nr(il) + 1 : cmap.get_nr(il) ;
    for ( int k=0; k<nk; ++k ) {
      set_vertex_ring( igl, cmap, xV, yV, fV ) ;
      ++igl ;
    }
  }
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF QUAD MESH
  //-----------------------------------------------------------------------------------------
  vector<int> regn_vlist, fR ;
  // ---
  int nR = 0 ;
  igl = 0 ;
  for ( int il=0; il<nl; ++il ) {
    // set number of rows of layer il
    int nr = cmap.get_nr(il) ;
    // il=0 : the first row is given by triangular cells
    if ( il==0 ) {
      set_tria_cells( igl, cmap, regn_vlist, fR, nR ) ;
      ++igl ;
      for ( int ic=1; ic<nr; ++ic ) {
	set_quad_cells( igl, cmap, regn_vlist, fR, nR ) ;
	++igl ;
      }
    } else {
      for ( int ic=0; ic<nr; ++ic ) {
	set_quad_cells( igl, cmap, regn_vlist, fR, nR ) ;
	++igl ;
      }
    }
  }
  // ---
  regn_vlist.push_back( nR ) ;
  //PRT(nR) ;
  //PRT( regn_vlist.size() ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR QUADRILATERAL-SHAPED MESH (underlying skeleton)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ;
  // set mesh name (standard flag 1xx)
  char str_flag[10] ; sprintf(str_flag,"%d",100+mesh_flag) ;
  char str_nr  [10] ; sprintf(str_nr,  "%d",cmap_flag)   ;
  char str_nt  [10] ; sprintf(str_nt,  "%d",nl)   ;
  string mesh_name = string("Quad-")+str_flag+string("-")+str_nr+string("x")+str_nt ;
  mesh . set_mesh_name( mesh_name ) ;

  //-----------------------------------------------------------------------------------------
  // set up curved face stuff
  //-----------------------------------------------------------------------------------------

#if CURVED_FACES
  mesh.setup_CoaxMap( coax ) ;
  mesh.setup_cvd_flist() ;

  bool ok_igl0(true), ok_igl1(true) ;

#define ALL_CVD_FACES 0

  int iR = 0 ;
  igl = 0 ;
  for ( int il=0; il<nl; ++il ) {
    // set number of rows of layer il
    int nr = cmap.get_nr(il) ;
    // il=0 : the first row is given by triangular cells
    if ( il==0 ) {
      set_cvd_tria_cells( igl, cmap, mesh, iR ) ;
      ++igl ;
      for ( int ic=1; ic<nr; ++ic ) {
#if ALL_CVD_FACES
	ok_igl0 = true ;
	ok_igl1 = true ;
#else
	// --- 
	//ok_igl0 = false ;    // curved faces at all internal layer interfaces
	//ok_igl1 = ic==nr-1 ; // and boundary
	// ---
	//ok_igl0 = false ;    // curved faces at all internal layer interfaces
	//ok_igl1 = ic==nr-1 && il!=nl-1 ; // NOT boundary
	// ---
	ok_igl0 = false ;      // curved faces only at the boundary
	ok_igl1 = ic==nr-1 && il==nl-1 ; 
#endif
	set_cvd_quad_cells( igl, cmap, mesh, iR, ok_igl0, ok_igl1 ) ;
	++igl ;
      }
    } else {
      for ( int ic=0; ic<nr; ++ic ) {
#if ALL_CVD_FACES
	ok_igl0 = true ;
	ok_igl1 = true ;
#else
	// ---
	//ok_igl0 = ic==0 ;    // curved faces at all internal layer interfaces
	//ok_igl1 = ic==nr-1 ; // and boundary
	// ---
	//ok_igl0 = ic==0 ;    // curved faces at all internal layer interfaces
	//ok_igl1 = ic==nr-1 && il!=nl-1 ; // NOT boundary
	// ---
	ok_igl0 = false ;      // curved faces only at the boundary
	ok_igl1 = ic==nr-1 && il==nl-1 ;
#endif
	set_cvd_quad_cells( igl, cmap, mesh, iR, ok_igl0, ok_igl1 ) ;
	++igl ;
      }
    }
  }
#endif

  //mesh.check_mesh() ;
}

// added ilev to the args list
// (original args list does not contain this input paramater)
void build_curved_face_mesh( mesh_2Dc & mesh, int nx=1, int ny=1, int mf=0,
			     double wx=0., double wy=0., int ilev=0 ) {

  int nl = ny ;     // used in case~1 as multiplying factor to refine the initial mesh
  switch ( mf/10 ) { 
  case 0 : // mesh flag 70x
    build_quad_cvd_mesh( mesh, nx, ny, mf, wx, wy, ilev ) ; break ; // 70x --> curved quads
    break ;
  case 1 : // mesh flag 71x
    build_disk_cvd_mesh( mesh, nx, nl, mf, wx, wy, ilev ) ; break ; // 71x --> curved sectors
    break ;
  default : assert(false) ;
  }
  // ---
}

#endif // end of _MESH_2D_CVD_BUILTIN_HH
