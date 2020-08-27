#ifndef _MESH_SUBGRID
#define _MESH_SUBGRID

// ==========================================================================================
/*
  
  SubGridCell contains the dataset to input the mesh builder

  sub_mesh2Dv_reader: reads input data from a file containing the sub-mesh description of the
                      subcells (base_mesh is not referred or used here) 
                      
  sub_mesh2Dv_painter: "paints" a line on a given base_mesh; output is stored in a vector of
                        SubGridCell objects

  sub_mesh2Dv_builder: build a sub-mesh

NOTA BENE NOTA BENE NOTA BENE NOTA BENE

  la versione di painter che usa MOF (solo per quads) e` nel file "subgrid_painter_MOF.hh"
  deve essere inglobata in questo file alla fine dello sviluppo

  la versione attuale di painter permette solo 2 materiali per cella, quindi le celle devono
  essere sufficientemente "piccole" per evitare inconsistenze (oppure si deve migliorare
  il codice)

 */
// ==========================================================================================

// simple structure to manage mesh input data
class SubGridCell {
public:
  int cell_id ;
  vector<int> regn_vlist, face_vlist, face_rlist ;
  vector<int> fV, fF, fR ;
  vector<double> xV, yV ;
public:
  SubGridCell( int _cell_id=0 ) : cell_id(_cell_id) {}
  ~SubGridCell() {}
} ;

// ===========================================================================================================================

class sub_mesh2Dv_reader {
private:
  int offset ;
  int n_subcells ;
  int cell_id, icell ;
  vector<SubGridCell> subcell_vec ;
  string fname ;

private:
  static istream & get_comments( istream & is ) ;
  void fatal_error( string fname ) ;
  void read_line   ( ifstream & inpf ) ;
  bool read_keywd  ( ifstream & inpf ) ;
  void read_subgrid( ifstream & inpf ) ;
  void read_subcell( ifstream & inpf ) ;
  void read_nodes  ( ifstream & inpf ) ;
  void read_regions( ifstream & inpf ) ;
  void read_faces  ( ifstream & inpf ) ;

public:
  sub_mesh2Dv_reader( string _fname ) : fname(_fname), icell(-1), offset(0) {}
  ~sub_mesh2Dv_reader() {}
  void read_file( string _fname="" ) ;
  void write_dataset() ;
  int n_subcell() { return subcell_vec.size() ; }
  SubGridCell & get_subcell( int ic ) { 
    assert( 0<=ic && ic<subcell_vec.size() ) ;
    return subcell_vec[ic] ; 
 }
} ;

istream & sub_mesh2Dv_reader::get_comments( istream & is ) {
  bool newline = true ;
  while( newline ) { 
    newline = false ;
    while( !is.eof() && (is.peek()==' ' || is.peek()=='\t' || is.peek()=='\n') ) { is.get() ; }
  }
  return is ;
}

void sub_mesh2Dv_reader::fatal_error( string fname ) {
  cerr << "fatal error in opening file " << fname << endl << flush ;
  assert(false) ;
}
void sub_mesh2Dv_reader::read_line( ifstream & inpf ) {
  string keywd = "" ;
  inpf >> get_comments >> keywd ;
  if ( keywd=="#" ) { 
    read_keywd( inpf ) ; 
  } else {
    //cout << keywd << endl ;
  }
}
bool sub_mesh2Dv_reader::read_keywd( ifstream & inpf ) {
  string keywd = "" ;
  inpf >> get_comments >> keywd ;
  // ----
  bool retval = true ;
  if      ( keywd == "subgrid" ) { read_subgrid(inpf) ; }
  else if ( keywd == "subcell" ) { read_subcell(inpf) ; }
  else if ( keywd == "nodes"   ) { read_nodes  (inpf) ; }
  else if ( keywd == "regions" ) { read_regions(inpf) ; }
  else if ( keywd == "faces"   ) { read_faces  (inpf) ; }
  else                           { retval = false     ; }
  return retval ;
}
void sub_mesh2Dv_reader::read_subgrid( ifstream & inpf ) {
  inpf >> n_subcells >> offset ;
  subcell_vec.resize(n_subcells) ;
  // ----
  //cout << "subgrid: n_subcells = " << n_subcells << endl ;
  //cout << endl ;
}
void sub_mesh2Dv_reader::read_subcell( ifstream & inpf ) {
  inpf >> cell_id ;
  ++icell ;
  subcell_vec[icell].cell_id = cell_id-offset ;
  // ----
  //cout << "subcell: cell_id = " << cell_id << endl ;
}
void sub_mesh2Dv_reader::read_nodes( ifstream & inpf ) {
  int nnodes, dum0, dum1, dum2, id ;
  double x, y ;
  inpf >> nnodes >> dum0 >> dum1 >> dum2 ;
  //cout << "subcell: nnodes = " << nnodes << endl ;
  for ( int i=0; i<nnodes; ++i ) {
    inpf >> id >> x >> y ;
    subcell_vec[icell].xV.push_back(x) ;
    subcell_vec[icell].yV.push_back(y) ;
    // -----
    //cout << id << "  " << x << "  " << y << endl ;
  }
}
void sub_mesh2Dv_reader::read_regions( ifstream & inpf ) {
  int nregions, nnodes, dum0, dum1, id, inode ;
  inpf >> nregions >> dum0 >> dum1 ;
  // ----
  //cout << "subcell: nregions = " << nregions << endl ;
  for ( int i=0; i<nregions; ++i ) {
    inpf >> id >> nnodes ;
    subcell_vec[icell].regn_vlist.push_back(nnodes) ;
    // -----
    //cout << id << "  " << nnodes << "  " ;
    for ( int j=0; j<nnodes; ++j ) {
      inpf >> inode ;
      subcell_vec[icell].regn_vlist.push_back(inode-offset) ;
      // -----
      //cout << inode << "  " ;
    }
    //cout << endl ;
  }
  subcell_vec[icell].regn_vlist.push_back(nregions) ;
}
void sub_mesh2Dv_reader::read_faces( ifstream & inpf ) {
  int nfaces, id, node0, node1, regn0, regn1, dum0 ;
  inpf >> nfaces >> dum0 ;
  // ----
  //cout << "subcell: nfaces = " << nfaces << endl ;
  for ( int i=0; i<nfaces; ++i ) {
    inpf >> id >> node0 >> node1 >> regn0 >> regn1 ;
    regn1 = regn1<0 ? UNSET+offset : regn1 ;
    subcell_vec[icell].face_vlist.push_back(node0-offset) ;
    subcell_vec[icell].face_vlist.push_back(node1-offset) ;
    subcell_vec[icell].face_rlist.push_back(regn0-offset) ;
    subcell_vec[icell].face_rlist.push_back(regn1-offset) ;
    // ----
    //cout << id << "  " << node0 << "  " << node1 << "  " << regn0 << "  " << regn1 << endl ;
  }
  //cout << endl ;
  subcell_vec[icell].face_vlist.push_back(nfaces) ;
  subcell_vec[icell].face_rlist.push_back(nfaces) ;
}
void sub_mesh2Dv_reader::read_file( string _fname ) {
  string inp_fname = _fname=="" ? fname : _fname ;
  ifstream inpf( inp_fname.c_str() ) ;
  if ( !inpf.good() ) { fatal_error(inp_fname) ; }
  while ( !inpf.eof() ) { read_line( inpf ) ; }
  inpf.close()  ;
}
void sub_mesh2Dv_reader::write_dataset() {
  string fname = string("subgrid.log") ;
  ofstream outp(fname.c_str()) ;
  outp << "subgrid: n_subcells = " << n_subcells << " offset = " << offset << endl ;
  for ( int ic=0; ic<n_subcells; ++ic ) {
    int cell_id = subcell_vec[ic].cell_id + offset ;
    outp << "subcell: cell_id = " << cell_id << endl ;
    // nodes
    int nnodes = subcell_vec[ic].xV.size() ;
    outp << "subcell: nnodes = " << nnodes << endl ;
    for ( int i=0; i<nnodes; ++i ) {
      int id = i+offset ;
      double x = subcell_vec[ic].xV[i] ;
      double y = subcell_vec[ic].yV[i] ;
      // -----
      outp << id << "  " << x << "  " << y << endl ;
    }
    // regions
    int nregions = subcell_vec[ic].regn_vlist.back() ;
    // ----
    outp << "subcell: nregions = " << nregions << endl ;
    int kk = 0 ;
    for ( int i=0; i<nregions; ++i ) {
      int id = i + offset ; 
      int nnodes = subcell_vec[ic].regn_vlist[kk++] ;
      // -----
      outp << id << "  " << nnodes << "  " ;
      for ( int j=0; j<nnodes; ++j ) {
	int inode = subcell_vec[ic].regn_vlist[kk++]+offset ;
	// -----
	outp << inode << "  " ;
      }
      outp << endl ;
    }
    // faces
    int nfaces = subcell_vec[ic].face_vlist.back() ;
    // ----
    outp << "subcell: nfaces = " << nfaces << endl ;
    int kv = 0 ;
    int kr = 0 ;
    for ( int i=0; i<nfaces; ++i ) {
      int id = i + offset ;
      int node0 = subcell_vec[ic].face_vlist[kv++]+offset ;
      int node1  =subcell_vec[ic].face_vlist[kv++]+offset ;
      int regn0 = subcell_vec[ic].face_rlist[kr++]+offset ;
      int regn1 = subcell_vec[ic].face_rlist[kr++]+offset ;
      // ----
      outp << id << "  " << node0 << "  " << node1 << "  " << regn0 << "  " << regn1 << endl ;
    }
    // final
    outp << endl ;
  }
  outp.close() ;
}

// ===========================================================================================================================

// "paint" a line on a mesh
class sub_mesh2Dv_painter {
private:
  static const int DIM = 2 ;
  static const int undef_zone = -1 ;
  const double epsi ;

  class face_isec {
  public:
    int iF ;
    double xF, yF, tF ;
    face_isec() : iF(UNSET), xF(0.), yF(0.), tF(0.) {}
    face_isec( int _iF, double _xF, double _yF, double _tF ) : iF(_iF), xF(_xF), yF(_yF), tF(_tF) {}
    ~face_isec() {}
  } ;

  vector<SubGridCell> subgridcell_vec ;

protected:
  double ax, ay, a0, tol ;
  mesh_2Dv & mesh ;

  vector<face_isec> face_isec_vec ;
  Vector face_isec_pnt ;
  VecInt face_isec_msk ;

  int  vrtx_izone( int iV ) ;
  void build_vrtx_izone_flag( Vector & vrtx_izone_vec ) ;

  void face_intersection( int iF, double & xF, double & yF, double & tF ) ;
  void face_intersection() ;

  void build_submesh() ;
  void build_submesh( int iR ) ;

  void setup() ;

public:
  sub_mesh2Dv_painter( mesh_2Dv & _mesh ) : mesh(_mesh), epsi(1e-12) {
    setup() ;
  }
  ~sub_mesh2Dv_painter() {}

  void paint_a_line( double _ax, double _ay, double _a0, double _tol ) ;
  void print_face_intersection() ;

  int n_subcell() { return subgridcell_vec.size() ; }
  SubGridCell & get_subcell( int ic ) {
    assert( 0<=ic && ic<subgridcell_vec.size() ) ;
    return subgridcell_vec[ic] ; 
  }
} ;

void sub_mesh2Dv_painter::setup() {
  int nV = mesh.n_vertex() ;
  int nF = mesh.n_face() ;
  face_isec_vec.resize(0) ;
  face_isec_pnt.setup(nF) ;
  face_isec_msk.setup(nF) ;
}

int sub_mesh2Dv_painter::vrtx_izone( int iV ) {
  assert( 0<=iV && iV<mesh.n_vertex() ) ;

  double xV = mesh.coords_V( iV, 0 ) ;
  double yV = mesh.coords_V( iV, 1 ) ;

  double rxy = ax*xV + ay*yV + a0 ;
  
  int retval = undef_zone ;         // undef  zone, probable intersection
  if ( rxy<-epsi ) { retval = 0 ; } // first  zone
  if ( rxy>+epsi ) { retval = 1 ; } // second zone
  
  return retval ;
}

void sub_mesh2Dv_painter::build_vrtx_izone_flag( Vector & vrtx_izone_vec ) {
  int nV = mesh.n_vertex() ;
  assert( vrtx_izone_vec.size()==nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    vrtx_izone_vec(iV) = vrtx_izone(iV) ;
  }
}

void sub_mesh2Dv_painter::face_intersection() {
  int nV = mesh.n_vertex() ;
  Vector vrtx_izone_vec(nV) ;
  build_vrtx_izone_flag( vrtx_izone_vec ) ;

  int isec = 0 ;
  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {
    
    int iV0 = mesh.face_vrtx(iF,0) ; 
    int iV1 = mesh.face_vrtx(iF,1) ; 

    int rV0 = vrtx_izone_vec(iV0) ;
    int rV1 = vrtx_izone_vec(iV1) ;

    bool check_intersection = rV0!=rV1 || rV0==undef_zone || rV1==undef_zone ;
    if ( check_intersection ) {
      double xF(0.), yF(0.), tF(0.) ;
      face_intersection( iF, xF, yF, tF ) ;
      if ( epsi<tF && tF<1.-epsi ) { // check internal intersection
	face_isec_vec.push_back( face_isec(iF,xF,yF,tF) ) ;
	face_isec_pnt(iF) = isec++ ;
	face_isec_msk(iF) = 1 ;
      }
    }
  }
}

void sub_mesh2Dv_painter::print_face_intersection() {
  int ns = face_isec_vec.size() ;
  for ( int is=0; is<ns; ++is ) {
    face_isec & fsec = face_isec_vec[is] ;
    int    iF = fsec.iF ;
    double xF = fsec.xF ;
    double yF = fsec.yF ;
    double tF = fsec.tF ;
    printf("iF = %3i xF = %14.7e  yF = %14.7e  tF=%14.7e\n",iF,xF,yF,tF) ;
  }
}

// r(x,y) = (ax,ay) * (x,y) + a0 = 0
// ( x(t), y(t) ) = v0 + t * ( v1-v0 )  
// tF = - ( a0 + (ax,ay)*(x,y) ) / ( (v1-v0)*(ax,ay) )  
// ( xF, yF ) = ( x(tF), y(tF) )
// DO NOT TREAT THE CASE WHERE THE FACE IS LOCATED ALONG THE INTERFACE
// THIS CASE IS REVEALED BY abs(den)>=1.e-12 AND THE ASSERT STOPS THE RUN
void sub_mesh2Dv_painter::face_intersection( int iF, double & xF, double & yF, double & tF ) {
  int    iV0 = mesh.face_vrtx(iF,0) ; 
  double xV0 = mesh.coords_V (iV0,0) ;
  double yV0 = mesh.coords_V (iV0,1) ;

  int    iV1 = mesh.face_vrtx(iF,1) ; 
  double xV1 = mesh.coords_V (iV1,0) ;
  double yV1 = mesh.coords_V (iV1,1) ;
  
  double den = ax*( xV1-xV0 ) + ay*( yV1-yV0 ) ;
  double num = ax*xV0 + ay*yV0 + a0 ;
  assert( abs(den)>1.e-12 ) ;

  // compute tF
  tF = -num/den ;
  if ( epsi<tF && tF<1.-epsi ) { // check internal intersection is non-degenerate
    if      ( tF<tol    ) { tF = tol ; }
    else if ( tF>1.-tol ) { tF = 1.-tol ; }
  }

  // compute xF, yF (intersection's coordinates)
  xF = xV0 + tF * ( xV1-xV0 ) ;
  yF = yV0 + tF * ( yV1-yV0 ) ;

  // debug printing
  //printf("qui: iF = %3i xF = %14.7e  yF = %14.7e  tF=%14.7e\n",iF,xF,yF,tF) ;
}

void sub_mesh2Dv_painter::build_submesh( int iR ) {

  int nRF = mesh.n_regn_face( iR ) ;
  int nRV = nRF ;

  int idx_V = 0 ;
  Vector idx_vrtx( mesh.n_vertex() ) ; // ??? sembra inutile....

  subgridcell_vec.push_back( SubGridCell(iR) ) ;
  SubGridCell & subc = subgridcell_vec.back() ;

  int new_nV = nRV+2 ;
  int new_nF = 0 ;
  vector<int> new_nodes ;

  int regn_count = 0 ;

  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = mesh.regn_face( iR, ilF ) ;
    
    int iV0 = mesh.ok_regn_face(iR,ilF) ? mesh.face_vrtx(iF,0) : mesh.face_vrtx(iF,1) ;
    int iV1 = mesh.ok_regn_face(iR,ilF) ? mesh.face_vrtx(iF,1) : mesh.face_vrtx(iF,0) ;

    subc.xV.push_back( mesh.coords_V(iV0,0) ) ;
    subc.yV.push_back( mesh.coords_V(iV0,1) ) ;
    subc.fV.push_back( -iV0-1 ) ;

    if ( face_isec_msk(iF)==1 ) {

      // additional node
      int ipsec = face_isec_pnt(iF) ;
      face_isec & fsec = face_isec_vec[ipsec] ;
      subc.xV.push_back( fsec.xF ) ;
      subc.yV.push_back( fsec.yF ) ;
      subc.fV.push_back( iF ) ;

      // first subface
      new_nF++ ;
      subc.face_vlist.push_back( idx_V   ) ;
      subc.face_vlist.push_back( (idx_V+1)%new_nV ) ;
      idx_vrtx(iV0) = idx_V++ ;
      // ---
      subc.face_rlist.push_back( regn_count ) ;
      subc.face_rlist.push_back( UNSET      ) ;
      // ---
      subc.fF.push_back(iF) ;

      // second subface
      new_nF++ ;
      new_nodes.push_back( idx_V ) ;
      subc.face_vlist.push_back( idx_V   ) ;
      subc.face_vlist.push_back( (idx_V+1)%new_nV ) ;
      idx_vrtx(iV0) = idx_V++ ;  // ??? sovrascrive il valore della prima faccia?
      // ---
      regn_count = (regn_count+1)%2 ;
      subc.face_rlist.push_back( regn_count ) ;
      subc.face_rlist.push_back( UNSET      ) ;
      // ---
      subc.fF.push_back(iF) ;

    } else {

      new_nF++ ;
      subc.face_vlist.push_back( idx_V   ) ;
      subc.face_vlist.push_back( (idx_V+1)%new_nV ) ;
      idx_vrtx(iV0) = idx_V++ ;
      // ---
      subc.face_rlist.push_back( regn_count ) ;
      subc.face_rlist.push_back( UNSET      ) ;
      // ---
      subc.fF.push_back(iF) ;
    }
  }

  // final face
  new_nF++ ;
  subc.face_vlist.push_back( new_nodes[0] ) ;
  subc.face_vlist.push_back( new_nodes[1] ) ;
  subc.face_vlist.push_back( new_nF ) ;
  // ----
  subc.face_rlist.push_back( 0 ) ;
  subc.face_rlist.push_back( 1 ) ;
  subc.face_rlist.push_back( new_nF ) ;
  // ----
  subc.fF.push_back( UNSET ) ;

#if 0 // useful for debugging
  int nV = idx_V ;
  for ( int iV=0; iV<nV; ++iV ) {
    printf("iV=%3i  xV[%2i]=%14.7e  yV[%2i]=%14.7e  fV[%2i]=%3i\n",iV,iV,subc.xV[iV],iV,subc.yV[iV],iV,subc.fV[iV]) ;
  }

  LINE(--) ;
  PRT_ARR(subc.face_vlist) ;

  LINE(--) ;
  PRT_ARR(subc.face_rlist) ;

  LINE(--) ;
  PRT_ARR(subc.fF) ;
#endif
}

void sub_mesh2Dv_painter::build_submesh() {
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    // check if the region is split in two subregions
    int nsec = 0 ;
    int nRF = mesh.n_regn_face( iR ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = mesh.regn_face( iR, ilF ) ;
      nsec += face_isec_msk(iF) ;
    }
    if ( nsec==2 ) { // found intersection
      build_submesh( iR ) ;
    }
  }
}

void sub_mesh2Dv_painter::paint_a_line( double _ax, double _ay, double _a0, double _tol ) {
  // ---
  ax  = _ax ;
  ay  = _ay ;
  a0  = _a0 ;
  tol = _tol ;
  // ---
  assert( ax>0. || abs(ax)<epsi && ay>0. ) ;
  // ---
  face_intersection() ; 
  build_submesh() ;
}

// ===========================================================================================================================

//#include "subgrid_painter_MOF.hh"
#include "subgrid_painter_MOF_nmats.hh"

// ===========================================================================================================================

class sub_mesh2Dv_builder {
  
private:
  mesh_2Dv & base_mesh ;
  
  // identify region and face origins
  void set_base_regn( int & cell_id, mesh_2Dv & submesh ) ;
  void set_base_face( int & cell_id, mesh_2Dv & submesh ) ;
  void renumber_mesh_faces( mesh_2Dv & submesh, vector<int> & face_perm ) ;
  
public:
  sub_mesh2Dv_builder( mesh_2Dv & _base_mesh ) : base_mesh(_base_mesh) {}
  ~sub_mesh2Dv_builder() {}

  void build_submesh( int base_iR, mesh_2Dv & submesh ) ;
  void build_submesh( int ic, sub_mesh2Dv_reader      & subgrid_reader,  mesh_2Dv & submesh ) ;
  void build_submesh( int ic, sub_mesh2Dv_painter_MOF & subgrid_painter, mesh_2Dv & submesh ) ;
  void build_submesh( int ic, sub_mesh2Dv_painter     & subgrid_painter, mesh_2Dv & submesh ) ;
} ;

// INTERNAL CONSTRUCTION of submesh (for cell base_iR)
void sub_mesh2Dv_builder::build_submesh( int base_iR, mesh_2Dv & submesh ) {

  // submesh parameter
#if 1
  int nx = 2 ;
  int ny = 3 ;
#else
  int nx = 2 + rand() % 3 ;
  int ny = 2 + rand() % 3 ;
#endif

  int mf = 0 ;
  double wx = 0 ;
  double wy = 0 ;

#if 0 // OLD
  // build the submesh on the [0,1]x[0,1] domain
  build_new_quad_mesh( submesh, nx, ny, mf, wx, wy ) ; // 1xx --> quads
  
  // get the vertex list of base_iR
  vector<int> regn_vlist ;
  base_mesh.get_regn_vrtx( base_iR, regn_vlist ) ;
  int nRV = regn_vlist.size() ;
  
  // reset the bounding box
  double x0(1e99), y0(1e99), x1(-1e99), y1(-1e99) ;
  for ( int ilV=0; ilV<nRV; ++ilV ) {
    int    iV = regn_vlist[ilV] ;
    double xV = base_mesh.coords_V( iV, 0 ) ;
    double yV = base_mesh.coords_V( iV, 1 ) ;
    x0 = min( x0, xV ) ;
    x1 = max( x1, xV ) ;
    y0 = min( y0, yV ) ;
    y1 = max( y1, yV ) ;
  }
  submesh.change_bbox( x0, y0, x1, y1 ) ;

#else // NEW, it works on any polygonal domain

  // build the submesh through a bilinear isoparametric transformation
  // of the [0,1]x[0,1] domain

  // get the vertex list of base_iR
  vector<int> regn_vlist ;
  base_mesh.get_regn_vrtx( base_iR, regn_vlist ) ;
  int nRV = regn_vlist.size() ;
  
  // reset the bounding box
  vector<double> xvrt, yvrt ;
  for ( int ilV=0; ilV<nRV; ++ilV ) {
    int    iV = regn_vlist[ilV] ;
    xvrt.push_back( base_mesh.coords_V( iV, 0 ) ) ;
    yvrt.push_back( base_mesh.coords_V( iV, 1 ) ) ;
  }

  // build the submesh on the [0,1]x[0,1] domain
  build_bil_quad_mesh( submesh, xvrt, yvrt, nx, ny, mf, wx, wy ) ; // 1xx --> quads

#endif

  // set index of base faces
  int cell_id = base_iR ;
  set_base_face( cell_id, submesh ) ;
}

// FROM EXTERNAL FILE
// this routine builds the submesh "ic"; 
// (i)    input data are from subgrid_reader; 
// (ii)   the initial mesh are without flags; then,
// (iii)  cell_id is determined by looking for the cell in 
//        base_mesh with the same barycenter
// (iV)   by looping on the boundary of submesh we assign 
//        the face ids of base_mesh corresponding to the boundary 
//        faces of submesh
void sub_mesh2Dv_builder::build_submesh( int ic, sub_mesh2Dv_reader & subgrid_reader, mesh_2Dv & submesh ) {

  // -------------------------------------
  // get subcell reference and region's id
  // -------------------------------------
  SubGridCell & subcell = subgrid_reader.get_subcell(ic) ;
  
  // -------------------------------------
  // work arrays
  // -------------------------------------
  vector<int>    fV, fF, face_vlist, face_rlist ;
  vector<double> xV, yV ;

  // -------------------------------------
  // set vertices
  // -------------------------------------
  int nV = subcell.xV.size() ;
  for ( int iV=0; iV<nV; ++iV ) {
    xV.push_back( subcell.xV[iV] ) ;
    yV.push_back( subcell.yV[iV] ) ;
    fV.push_back( UNSET ) ;
  }
  
  // -------------------------------------
  // set faces
  // -------------------------------------
  int nl = subcell.face_vlist.size() ;
  for ( int il=0; il<nl; ++il ) {
    face_vlist.push_back( subcell.face_vlist[il] ) ;
    face_rlist.push_back( subcell.face_rlist[il] ) ;
  }
  for ( int il=0; il<face_vlist.back(); ++il ) {
    fF.push_back( UNSET ) ;
  }

  // -------------------------------------
  // build the mesh
  // -------------------------------------
  mesh2Dv_builder mesh_builder(submesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, face_vlist, face_rlist, fF ) ;

  // -------------------------------------
  // REQUIRED: identify the base cell 
  // -------------------------------------
  int cell_id = subcell.cell_id ;
  set_base_regn( cell_id, submesh ) ;
  set_base_face( cell_id, submesh ) ;
}

// FROM PAINTER (interface input form data.lap)
void sub_mesh2Dv_builder::build_submesh( int ic, sub_mesh2Dv_painter & subgrid_painter, mesh_2Dv & submesh ) {

  // -------------------------------------
  // get subcell reference and region's id
  // -------------------------------------
  SubGridCell & subcell = subgrid_painter.get_subcell(ic) ;
  
  // -------------------------------------
  // work arrays
  // -------------------------------------
  vector<int>    fV, fF, face_vlist, face_rlist ;
  vector<double> xV, yV ;

  // -------------------------------------
  // set vertices
  // -------------------------------------
  int nV = subcell.xV.size() ;
  for ( int iV=0; iV<nV; ++iV ) {
    xV.push_back( subcell.xV[iV] ) ;
    yV.push_back( subcell.yV[iV] ) ;
    fV.push_back( UNSET ) ;
  }
  
  // -------------------------------------
  // set faces
  // -------------------------------------
  int nl = subcell.face_vlist.size() ;
  for ( int il=0; il<nl; ++il ) {
    face_vlist.push_back( subcell.face_vlist[il] ) ;
    face_rlist.push_back( subcell.face_rlist[il] ) ;
  }
  for ( int il=0; il<face_vlist.back(); ++il ) {
    fF.push_back( subcell.fF[il] ) ;
  }

  // -------------------------------------
  // build the mesh
  // -------------------------------------
  mesh2Dv_builder mesh_builder(submesh) ;
  mesh_builder . build_the_mesh( xV, yV, fV, face_vlist, face_rlist, fF ) ;

  // -------------------------------------
  // set the mesh name
  // -------------------------------------
  string str_cell_id("subgrid: cell_id ") ;
  str_cell_id = append_int_to_string( subcell.cell_id, str_cell_id ) ;
  submesh.set_mesh_name( str_cell_id ) ;

  // -------------------------------------
  // check cell_id 
  // (base face indices have been 
  // already provided by painter)
  // -------------------------------------
  int cell_id = subcell.cell_id ;
  set_base_regn( cell_id, submesh ) ;
  //set_base_face( cell_id, submesh ) ;
}

// FROM PAINTER_MOF (interface input form data.lap)
void sub_mesh2Dv_builder::build_submesh( int ic, sub_mesh2Dv_painter_MOF & subgrid_painter, mesh_2Dv & submesh ) {
  MSGF("begin  sub_mesh2Dv_builder::build_submesh") ;

  // -------------------------------------
  // get subcell reference and region's id
  // -------------------------------------
  SubGridCell & subcell = subgrid_painter.get_subcell(ic) ;

  // -------------------------------------
  // build the mesh
  // -------------------------------------
  mesh2Dv_builder mesh_builder(submesh) ;
  mesh_builder . build_the_mesh( subcell.xV, subcell.yV, subcell.fV, subcell.regn_vlist, subcell.fR ) ;

  // -------------------------------------
  // set the mesh name
  // -------------------------------------
  string str_cell_id("subgrid: cell_id ") ;
  str_cell_id = append_int_to_string( subcell.cell_id, str_cell_id ) ;
  submesh.set_mesh_name( str_cell_id ) ;

  // -------------------------------------
  // renumber the mesh faces (boundary faces first)
  // -------------------------------------
  vector<int> face_perm ;
  renumber_mesh_faces( submesh, face_perm ) ;
  mesh_builder.reorder_faces( face_perm ) ;

  // -------------------------------------
  // check cell_id 
  // (base face indices have been 
  // already provided by painter)
  // -------------------------------------
  int cell_id = subcell.cell_id ;
  set_base_regn( cell_id, submesh ) ;
  set_base_face( cell_id, submesh ) ;

  MSGF("end of sub_mesh2Dv_builder::build_submesh") ;
}

// aux method: determine/check cell_id
void sub_mesh2Dv_builder::set_base_regn( int & cell_id, mesh_2Dv & submesh ) {
  // compute the barycenter of the submesh
  int nR = submesh.n_region() ;
  double xR(0.), yR(0.), mR(0) ;
  for ( int iR=0; iR<nR; ++iR ) {
    double aR = submesh.get_regn_measure(iR) ; 
    xR += aR * submesh.coords_R(iR,0) ;
    yR += aR * submesh.coords_R(iR,1) ;
    mR += aR ;
  }
  xR /= mR ;
  yR /= mR ;

  // look for the cell in the base mesh
  int base_cell_id = UNSET ;
  double tol = 1.e-6 * sqrt( xR*xR+yR*yR ) ;
  int base_nR = base_mesh.n_region() ;
  bool do_loop(true) ;
  for ( int base_iR=0; base_iR<base_nR && do_loop ; ++base_iR ) {
    double diff_R = sqrt( pow(xR-base_mesh.coords_R(base_iR,0),2) + 
			  pow(yR-base_mesh.coords_R(base_iR,1),2) ) ;
    if ( diff_R<tol ) {
      do_loop = false ;
      base_cell_id = base_iR ;
    }
  }

  // check a cell_id was found
  assert( !do_loop ) ;

  // check consistency if an input cell_id is given
  if ( cell_id!=UNSET ) {
    if ( !base_cell_id==cell_id ) {
      cout << "cell_id from input file is not consistent with base_mesh " << endl ;
      cout << "cell_id (from base_mesh)  = " << base_cell_id << endl ;
      cout << "cell_id (from input data) = " << cell_id << endl ;
      assert( base_cell_id==cell_id ) ;
    }
  }
}

// set index face from base mesh to faces of submesh
void sub_mesh2Dv_builder::set_base_face( int & cell_id, mesh_2Dv & submesh ) {
  // look for the boundary edges (the abs of pscal is required because
  // to take care of the mutual orientation of the normal vectors)
  int base_iR = cell_id ;

  double tol = 1.e-6 * sqrt( pow(base_mesh.coords_R(cell_id,0),2) + 
			     pow(base_mesh.coords_R(cell_id,1),2) ) ;
  int nbF = submesh.n_bface() ;
  int base_nRF = base_mesh.n_regn_face( base_iR ) ;
  for ( int ilF=0; ilF<nbF; ++ilF ) {
    int    iF  = submesh.get_bnd_face(ilF) ;
    double nxF = submesh.get_nor(iF,0) ;
    double nyF = submesh.get_nor(iF,1) ;
    double xF  = submesh.coords_F(iF,0) ;
    double yF  = submesh.coords_F(iF,1) ;
    bool do_loop(true) ;
    for ( int klF=0; klF<base_nRF && do_loop; ++klF ) {
      int base_iF  = base_mesh.regn_face(base_iR,klF) ;
      double base_nxF = base_mesh.get_nor(base_iF,0) ;
      double base_nyF = base_mesh.get_nor(base_iF,1) ;
      double pscal = abs(nxF*base_nxF+nyF*base_nyF) ;
      if ( abs(1.-pscal)<tol ) {
	double base_xF = base_mesh.coords_F(base_iF,0) ;                    // the cross product is to deal
	double base_yF = base_mesh.coords_F(base_iF,1) ;                    // for the case when the base face 
	double check_F = abs( (xF-base_xF)*base_mesh.get_tng(base_iF,1) -   // is split into two subfaces
			      (yF-base_yF)*base_mesh.get_tng(base_iF,0) ) ; // 
	if ( check_F<tol ) {                                                    
	  do_loop = false ;
	  submesh.set_fF(iF,base_iF) ;
	}
      }
    } // end of --->> for ( int klF=0; klF<base_nRF && do_loop; ++klF ) {...
    assert( !do_loop ) ;
  } // end of --->> for ( int ilF=0; ilF<nbF; ++ilF ) {...
}

void sub_mesh2Dv_builder::renumber_mesh_faces( mesh_2Dv & submesh, vector<int> & face_perm ) {
  // ----
  int nF = submesh.n_face() ;  
  face_perm.resize(nF) ;
  // ----
  int new_iF(0) ;  
  for ( int iF=0; iF<nF; ++iF ) {
    if ( submesh.is_boundary_face(iF) ) {
      face_perm[iF] = new_iF++ ;
    }
  }
  for ( int iF=0; iF<nF; ++iF ) {
    if ( submesh.is_internal_face(iF) ) {
      face_perm[iF] = new_iF++ ;
    }
  }
}

// ===========================================================================================================================

class submesh_post_proc {
private:
  mesh_2Dv          & base_mesh ;
  vector<mesh_2Dv*> & submesh_vec ;

public:
  submesh_post_proc( mesh_2Dv & _base_mesh, vector<mesh_2Dv*> & _submesh_vec ) : 
    base_mesh(_base_mesh), submesh_vec(_submesh_vec) {}
  ~submesh_post_proc() {}
  
  void plot_all( string mesh_name, string psout_file, bool draw_label=false ) ;

} ;

void submesh_post_proc::plot_all( string mesh_name, string psout_file, bool draw_label ) {
  MSGF("begin  submesh_post_proc::plot_all") ;
  
  PS_Driver psd(psout_file) ;
  mesh2D_psout msh_psout( base_mesh, psd ) ;

  // just to change offset (if required)
  msh_psout.set_offset( OFFSET ) ;

  // draw the base mesh
  bool draw_numbers = false ;
  msh_psout . draw_mesh(string("mesh"),draw_numbers) ;
  
  // draw the sub-grid meshes
  int n_submesh = submesh_vec.size() ;
  for ( int isub=0; isub<n_submesh; ++isub ) {
    mesh_2Dv & submesh = *(submesh_vec[isub]) ;
    string str_cell_id = string("% ") + submesh.get_mesh_name() ;
    str_cell_id = append_int_to_string( isub, str_cell_id ) ;
    psd.dump_string(str_cell_id) ;
    psd.set_rgb_colour( "MAGENTA" ) ;
    //psd.set_line_type ( "DASH1" ) ;
    psd.set_line_type ( "SOLID" ) ;
    for ( int iF=0; iF<submesh.n_face(); ++iF ) {
      if ( submesh.is_internal_face(iF) ) {
	// draw the submesh
	int    iV0   = submesh.face_vrtx(iF,0) ;
	int    iV1   = submesh.face_vrtx(iF,1) ;
	double xV[2] = { submesh.coords_V(iV0,0), submesh.coords_V(iV1,0) } ;
	double yV[2] = { submesh.coords_V(iV0,1), submesh.coords_V(iV1,1) } ;
	psd.draw_line( xV[0], yV[0], xV[1], yV[1] ) ;
      }
    }
  }
  
  if ( draw_label ) {
    for ( int isub=0; isub<n_submesh; ++isub ) {
      //int cell_id = ...
      mesh_2Dv & submesh = *(submesh_vec[isub]) ;
      
      // draw vertex numbers
      psd.dump_string( "% ----- vrtx numbers" ) ;
      psd.set_rgb_colour( "BLACK" ) ;
      psd.set_font( "Vrtx" ) ;
      for ( int iV=0; iV<submesh.n_vertex(); ++iV ) {
	double xV = submesh.coords_V(iV,0) ;
	double yV = submesh.coords_V(iV,1) ;
	psd.draw_number( iV+msh_psout.get_offset(), xV, yV ) ;
      }
      
      for ( int iF=0; iF<submesh.n_face(); ++iF ) {
	// draw edge numbers
	psd.dump_string( "% ----- edge numbers" ) ;
	psd.set_rgb_colour( "RED" ) ;
	psd.set_font( "Edge" ) ;
	double xF = submesh.coords_F(iF,0) ;
	double yF = submesh.coords_F(iF,1) ;
	psd.draw_number( iF+msh_psout.get_offset(), xF, yF ) ;
      }
      
      // draw regn numbers
      psd.dump_string( "% ----- regn numbers" ) ;
      psd.set_rgb_colour( "BLUE" ) ;
      psd.set_font( "Poly" ) ;
      for ( int iR=0; iR<submesh.n_region(); ++iR ) {
	double xR = submesh.coords_R(iR,0) ;
	double yR = submesh.coords_R(iR,1) ;
	psd.draw_number( iR+msh_psout.get_offset(), xR, yR ) ;
      }
    }
  }

  MSGF("end of submesh_post_proc::plot_all") ;
}

#endif // end of _MESH_SUBGRID
