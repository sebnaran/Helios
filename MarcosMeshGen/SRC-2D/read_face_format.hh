/* example:
## the first line is always a comment!

# offset 0

# nodes
6 2 0 0
   0   0.0   0.0
   1   0.5   0.0
   2   1.0   0.0
   3   0.0   1.0
   4   0.5   1.0
   5   1.0   1.0

# regions
2 0 0
   0   4    0  1  4  3 
   1   4    1  2  5  4

# faces
7 0 
   0   0  1   0  -1
   1   3  0   0  -1
   2   1  2   1  -1
   3   1  4   0   1
   4   2  5   1  -1
   5   4  3   0  -1
   6   5  4   1  -1
*/

/*
  
  Reading of flags is not yet implemented. Internal flags are all set to UNSET.

 */

class mesh2D_reader_FaceFormat : public mesh2D_reader {

private:
  static istream & get_comments( istream & is ) {
    bool newline = true ;
    while( newline ) { 
      newline = false ;
      while( !is.eof() && (is.peek()==' ' || is.peek()=='\t' || is.peek()=='\n') ) { is.get() ; }
    }
    return is ;
  }  

  void fatal_error( string fname ) {
    cerr << "fatal error in opening file " << fname << endl << flush ;
    assert(false) ;
  }

private:
  mesh_2Dv & inp_mesh  ;
  string     file_name ;
  int        offset    ;    // ==1 fortran offset

  // instantiate work arrays
  vector<int>    face_vlist, face_rlist, regn_vlist, fV, fF, fR ; 
  vector<double> xV, yV ; 
  
  // methods
  void read_nodes  ( ifstream & inpf ) ;
  void read_regions( ifstream & inpf ) ;
  void read_faces  ( ifstream & inpf ) ;
  void read_file   () ;
  void read_line   ( ifstream & inpf ) ;
  void read_offset ( ifstream & inpf ) ;
  bool read_keywd  ( ifstream & inpf ) ;

  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV ) { assert(false) ; }
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) { assert(false) ; }
  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR ) {
    read_file() ;
  }

public:
  mesh2D_reader_FaceFormat( mesh_2Dv & _inp_mesh, string _file_name, int _offset=OFFSET ) : 
    inp_mesh(_inp_mesh), file_name(_file_name), offset(_offset) {}
  ~mesh2D_reader_FaceFormat() {}

  void read_and_build() ;
} ;
void mesh2D_reader_FaceFormat::read_and_build() {
  DBGF(begin  -->> mesh2D_reader_FaceFormat) ;
  
  // input mesh from data file
  read_file() ;

  for ( int i=0; i<xV.size(); ++i ) {
    printf("i=%i  xV[%i]=%14.7e  yV[%i]=%14.7e  fV[%i]=%i\n",i,i,xV[i],i,yV[i],i,fV[i]) ;
  }

  for ( int i=0; i<face_vlist.size(); ++i ) {
    printf("i=%i  face_vlist[%i]=%2i  face_rlist[%i]=%2i\n",i,i,face_vlist[i],i,face_rlist[i]) ;
  }
  

  // start builder
  mesh2Dv_builder mesh_builder( inp_mesh ) ;
  mesh_builder . build_the_mesh( xV,yV, fV, face_vlist, face_rlist, fF ) ;

  DBGF(end of -->> mesh2D_reader_FaceFormat) ;
}
void mesh2D_reader_FaceFormat::read_nodes( ifstream & inpf ) {
  int nnodes, dum0, dum1, dum2, id ;
  double x, y ;
  inpf >> nnodes >> dum0 >> dum1 >> dum2 ;
  cout << "nnodes = " << nnodes << endl ;
  for ( int i=0; i<nnodes; ++i ) {
    inpf >> id >> x >> y ;
    xV.push_back(x) ;
    yV.push_back(y) ;
    fV.push_back(UNSET) ;
    // -----
    cout << id << "  " << x << "  " << y << endl ;
  }
}
void mesh2D_reader_FaceFormat::read_regions( ifstream & inpf ) {
  int nregions, nnodes, dum0, dum1, id, inode ;
  inpf >> nregions >> dum0 >> dum1 ;
  // ----
  cout << "nregions = " << nregions << endl ;
  for ( int i=0; i<nregions; ++i ) {
    inpf >> id >> nnodes ;
    regn_vlist.push_back(nnodes) ;
    // -----
    cout << id << "  " << nnodes << "  " ;
    for ( int j=0; j<nnodes; ++j ) {
      inpf >> inode ;
      regn_vlist.push_back(inode-offset) ;
      // -----
      cout << inode << "  " ;
    }
    cout << endl ;
    fR.push_back(UNSET) ;
  }
  regn_vlist.push_back(nregions) ;
}
void mesh2D_reader_FaceFormat::read_faces( ifstream & inpf ) {
  int nfaces, id, node0, node1, regn0, regn1, dum0 ;
  inpf >> nfaces >> dum0 ;
  // ----
  cout << "nfaces = " << nfaces << endl ;
  for ( int i=0; i<nfaces; ++i ) {
    inpf >> id >> node0 >> node1 >> regn0 >> regn1 ;
    regn1 = regn1<0 ? UNSET+offset : regn1 ;
    face_vlist.push_back(node0-offset) ;
    face_vlist.push_back(node1-offset) ;
    face_rlist.push_back(regn0-offset) ;
    face_rlist.push_back(regn1-offset) ;
    // ----
    fF.push_back(UNSET) ;
    cout << id << "  " << node0 << "  " << node1 << "  " << regn0 << "  " << regn1 << endl ;
  }
  //cout << endl ;
  face_vlist.push_back(nfaces) ;
  face_rlist.push_back(nfaces) ;
}
void mesh2D_reader_FaceFormat::read_file() {
  ifstream inpf( file_name.c_str() ) ;
  if ( !inpf.good() ) { fatal_error(file_name) ; }
  while ( !inpf.eof() ) { read_line( inpf ) ; }
  inpf.close()  ;
}
void mesh2D_reader_FaceFormat::read_line( ifstream & inpf ) {
  string keywd = "" ;
  inpf >> eatcomments >> keywd ;
  PRT(keywd) ;
  // ----
  if ( keywd=="#" ) { 
    read_keywd( inpf ) ; 
  } else {
    //cout << keywd << endl ;
  }
}
void mesh2D_reader_FaceFormat::read_offset( ifstream & inpf ) {
  inpf >> offset ;
  // ----
  cout << "offset = " << offset << endl ;
  cout << endl ;
}
bool mesh2D_reader_FaceFormat::read_keywd( ifstream & inpf ) {
  string keywd = "" ;
  inpf >> get_comments >> keywd ;
  // ----
  bool retval = true ;
  if      ( keywd == "offset"  ) { read_offset (inpf) ; }
  else if ( keywd == "nodes"   ) { read_nodes  (inpf) ; }
  else if ( keywd == "regions" ) { read_regions(inpf) ; }
  else if ( keywd == "faces"   ) { read_faces  (inpf) ; }
  else                           { retval = false     ; }
  return retval ;
}
