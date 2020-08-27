#ifndef _MESH_2D_GREADER_HH
#define _MESH_2D_GREADER_HH

#include <iostream>
#include<fstream>
#include <vector>
#include <cassert>
#include<string>
using namespace std ;

void print_logfile( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
		    vector<int> & regn_vlist, vector<int> & fR, int offset=OFFSET ) {
  
  // open log file
  ofstream LOGF("mesh.log") ;
  
  // print vertex data
  int nV = fV.size() ;
  LOGF << "number of vertices " << nV << endl ;
  for ( int iV=0; iV<nV; ++iV ) {
    LOGF << iV+offset << "  " << xV[iV] << "  " << yV[iV] << "  " << "  " << fV[iV] << endl ;
  }
  LOGF << "---------------------------------" << endl ;
  
  // print regn data
  int nR = regn_vlist.back() ;
  LOGF << "number of regions " << nR << endl ;
  int kV = 0;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRV = regn_vlist[kV++] ;
    LOGF << nRV << "\t< " ;
    for ( int iRV=0; iRV<nRV-1; ++iRV ) {
      LOGF << regn_vlist[kV++]+offset << ", " ;
    }
    LOGF << regn_vlist[kV++]+offset << " >  " ;
    LOGF << fR[iR] << " " << endl ;
  }

  // close log file
  LOGF.close() ;
}

//--------------------------------------------------------------------------------------------
// base class for 2D mesh readers to be used in public derivations
//--------------------------------------------------------------------------------------------

class mesh2D_reader {

protected:
  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }
  
  static istream & eatchar(istream & s) { s.get() ; return s ; }
  
  static istream & eatcomments(istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
	    && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }
  
  void error_message( string & s ) {
    cerr << "fatal error:\n" 
	 << "mesh_reader --> read_mesh() cannot open file " << s << endl ;
    exit(0) ;
  }

  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV )=0 ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR )=0 ;
  
public:
  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR )=0 ;
} ;

//--------------------------------------------------------------------------------------------
// read input of MESH_2D in "General Format", which is a generalization of TRIANGLE format
//--------------------------------------------------------------------------------------------

class mesh2D_reader_GeneralFormat : public mesh2D_reader {
  
private:
  mesh_2Dv & inp_mesh  ;
  string     file_name ;
  int        offset    ;    // ==1 fortran offset

  int nodelem ;
  void read_regn_gfree_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_gfree_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_fixed_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_fixed_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  
protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) ;
  
public:
  mesh2D_reader_GeneralFormat( mesh_2Dv & _inp_mesh, string _file_name, int _offset=OFFSET ) : 
    inp_mesh(_inp_mesh), file_name(_file_name), offset(_offset) {}
  ~mesh2D_reader_GeneralFormat() {}

  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR ) ;
  
  void read_and_build() ;
} ;
void mesh2D_reader_GeneralFormat :: read_and_build() {
  DBGF(begin  -->> mesh2D_reader_GeneralFormat) ;
  

  // instantiate work arrays
  vector<int>    fV, regn_vlist, fR ; 
  vector<double> xV, yV ; 
  
  // input mesh from data file
  read_the_mesh( xV, yV, fV, regn_vlist, fR         ) ;
  //print_logfile( xV, yV, fV, regn_vlist, fR, offset ) ;

  // start builder
  mesh2Dv_builder mesh_builder( inp_mesh ) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ; // !!!
  
  DBGF(end of -->> mesh2D_reader_GeneralFormat) ;
}
void mesh2D_reader_GeneralFormat ::read_the_mesh( vector<double> & xV, 
						  vector<double> & yV,  
						  vector<int>    & fV, 
						  vector<int>    & regn_vlist, 
						  vector<int>    & fR ) {
  // read NODE file of triangle output
  string node_file_name = file_name + ".node" ;
  ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;
  
  // read ELE file of triangle output
  string ele_file_name = file_name + ".ele" ;
  ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_vlist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
}
void mesh2D_reader_GeneralFormat :: read_vrtx_data ( ifstream & input_file, 
						     vector<double> & xV, 
						     vector<double> & yV,
						     vector<int>    & fV ) {
  MSG("start mesh2D_reader_GeneralFormat::read_vrtx_data"<<endl) ;
  
  // declarations
  double x(0.), y(0.) ;
  int id(0), nV(0), ndim(0), n_attr(0), nbmrk(0), f(0) ;
  
  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  
  // resize the vertex list of mesh
  MSG("---#vertices: ") ; PRT( nV ) ;
  // read mesh vertices
  if ( nbmrk==0 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( UNSET ) ;
    }
  } else if ( nbmrk==1 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( f ) ;
    } 
  } else {
    assert(0) ;
  }
  MSG("end of mesh_reader_GeneralFormat::read_vrtx_data"<<endl) ;
}
void mesh2D_reader_GeneralFormat :: 
read_regn_gfree_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), nRV(0), kV(0) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( int k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( UNSET ) ;
  }
}
void mesh2D_reader_GeneralFormat ::
read_regn_gfree_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), nRV(0), kV(0), regn_flag(UNSET) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( int k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat :: 
read_regn_fixed_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), kV(0) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( int k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( UNSET ) ;
  }
}
void mesh2D_reader_GeneralFormat :: 
read_regn_fixed_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), kV(0), regn_flag(UNSET) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( int k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat :: read_regn_data( ifstream    & ele_file, 
						    vector<int> & regn_vlist,
						    vector<int> & fR ) {
  MSG("start read_regn_data"<<endl) ;

  // declarations
  int nR(0), nbmrk(0) ;
  
  // read file header
  ele_file >> eatcomments >> nR >> nodelem >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  MSG("---#regions: ") ; PRT(nR) ;
  if      ( nodelem==0 && nbmrk==0 ) { read_regn_gfree_format_int_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem==0 && nbmrk==1 ) { read_regn_gfree_format_ext_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==0 ) { read_regn_fixed_format_int_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==1 ) { read_regn_fixed_format_ext_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else { assert(0) ; }
  regn_vlist.push_back( nR ) ;

  MSG("end of read_regn_data"<<endl) ;
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class mesh2D_reader_QVoronoiOutput : public mesh2D_reader {
  
private:
  mesh_2Dv & inp_mesh  ;
  string     file_name ;
  int        offset    ;    // ==1 fortran offset
  int        nDim, nV, nR, nmR ; // nmR = number meaningful regions
  
protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) ;
  
public:
  mesh2D_reader_QVoronoiOutput( mesh_2Dv & _inp_mesh, string _file_name, int _offset=1 ) : 
    inp_mesh(_inp_mesh), file_name(_file_name), offset(1) {}
  ~mesh2D_reader_QVoronoiOutput() {}

  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR ) ;
  
  void read_and_build() ;
} ;
void mesh2D_reader_QVoronoiOutput :: read_and_build() {
  DBGF(begin  -->> mesh2D_reader_QVoronoiOutput) ;
  

  // instantiate work arrays
  vector<int>    fV, regn_vlist, fR ; 
  vector<double> xV, yV ; 
  
  // input mesh from data file
  read_the_mesh( xV, yV, fV, regn_vlist, fR         ) ;
  print_logfile( xV, yV, fV, regn_vlist, fR, offset ) ;

  for ( int i=0; i<xV.size(); ++i ) {
    VAL(i) ; VALA(i,xV) ; VALA(i,yV); PRTA(i,fV) ;
  }
  LINE(--) ;
  PRT_ARR(regn_vlist) ; 
  LINE(--) ;
  PRT_ARR(fR) ; 
  
  // start builder
  mesh2Dv_builder mesh_builder( inp_mesh ) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ; // !!!
  
  DBGF(end of -->> mesh2D_reader_QVoronoiOutput) ;
}
void mesh2D_reader_QVoronoiOutput ::read_the_mesh( vector<double> & xV, 
						   vector<double> & yV,  
						   vector<int>    & fV, 
						   vector<int>    & regn_vlist, 
						   vector<int>    & fR ) {
  // read DATA file
  string data_file_name = file_name + ".data" ;
  ifstream data_file( data_file_name.c_str() ) ;
  if ( data_file.good() ) {
    // read the header
    data_file >> nDim >> eatline >> eatcomments >> nV >> nR >> nmR >> eatline ;
    VAL(nDim) ; VAL(nV) ; VAL(nR) ; PRT(nmR) ;
    // read the vertex data
    read_vrtx_data( data_file, xV, yV, fV ) ;
    // read the region data
    read_regn_data( data_file, regn_vlist, fR ) ;
  } else {
    error_message(data_file_name) ;
  }
  data_file.close() ;
}
void mesh2D_reader_QVoronoiOutput :: read_vrtx_data ( ifstream & input_file, 
						      vector<double> & xV, 
						      vector<double> & yV,
						      vector<int>    & fV ) {
  MSG("begin mesh2D_reader_QVoronoiOutput::read_vrtx_data"<<endl) ;
  
  // declarations
  double x(0.), y(0.) ;
  
  // skip the infinity point 
  input_file >> eatline ;
  
  // resize the vertex list of mesh
  MSG("---#vertices: ") ; PRT( nV ) ;
  for ( int iV=1 ; iV<nV ; ++iV ) {
    input_file >> eatcomments >> x >> y >> eatline ;
    xV.push_back( x ) ;
    yV.push_back( y ) ;
    fV.push_back( UNSET ) ;
  }
  MSG("end of mesh_reader_QVoronoiOutput::read_vrtx_data"<<endl) ;
}
void mesh2D_reader_QVoronoiOutput :: read_regn_data( ifstream & input_file, 
						     vector<int> & regn_vlist, 
						     vector<int> & fR ) {
  MSG("begin mesh2D_reader_QVoronoiOutput::read_regn_data"<<endl) ;
  int nRV(0), kV(0) ;
  PRT(nR) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    input_file >> eatcomments >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( int k=0; k<nRV; ++k ) {
      input_file >> kV ;
      regn_vlist.push_back( kV - offset) ; // offset = 1 as we skip the infinity node
    }
    input_file >> eatline ;
    fR.push_back( UNSET ) ;
  }
  regn_vlist.push_back(nR) ;
  MSG("end of mesh2D_reader_QVoronoiOutput::read_regn_data"<<endl) ;
}

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
  Notes: 
  
  1) Reading of flags is not yet implemented. Internal flags are all
  set to UNSET.

  2) Regions are input, but actually not used to generate the mesh. The
  block "regions" can be omitted in the input file.

  3) This implementation is very similar to the one given in
  mesh2D_subgrid.hh to input the grid associated with subcell
  substructures.

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
  //DBGF(begin  -->> mesh2D_reader_FaceFormat) ;
  
  // input mesh from data file
  read_file() ;

#if 0 // debugging, TBR
  for ( int i=0; i<xV.size(); ++i ) {
    printf("i=%i  xV[%i]=%14.7e  yV[%i]=%14.7e  fV[%i]=%i\n",i,i,xV[i],i,yV[i],i,fV[i]) ;
  }
  for ( int i=0; i<face_vlist.size(); ++i ) {
    printf("i=%i  face_vlist[%i]=%2i  face_rlist[%i]=%2i\n",i,i,face_vlist[i],i,face_rlist[i]) ;
  }
#endif

  // start builder
  mesh2Dv_builder mesh_builder( inp_mesh ) ;
  mesh_builder . build_the_mesh( xV,yV, fV, face_vlist, face_rlist, fF ) ;

  //DBGF(end of -->> mesh2D_reader_FaceFormat) ;
}
void mesh2D_reader_FaceFormat::read_nodes( ifstream & inpf ) {
  int nnodes, dum0, dum1, dum2, id ;
  double x, y ;
  inpf >> nnodes >> dum0 >> dum1 >> dum2 ;
  //cout << "nnodes = " << nnodes << endl ;
  for ( int i=0; i<nnodes; ++i ) {
    inpf >> id >> x >> y ;
    xV.push_back(x) ;
    yV.push_back(y) ;
    fV.push_back(UNSET) ;
    // -----
    //cout << id << "  " << x << "  " << y << endl ;
  }
}
void mesh2D_reader_FaceFormat::read_regions( ifstream & inpf ) {
  int nregions, nnodes, dum0, dum1, id, inode ;
  inpf >> nregions >> dum0 >> dum1 ;
  // ----
  //cout << "nregions = " << nregions << endl ;
  for ( int i=0; i<nregions; ++i ) {
    inpf >> id >> nnodes ;
    regn_vlist.push_back(nnodes) ;
    // -----
    //cout << id << "  " << nnodes << "  " ;
    for ( int j=0; j<nnodes; ++j ) {
      inpf >> inode ;
      regn_vlist.push_back(inode-offset) ;
      // -----
      //cout << inode << "  " ;
    }
    //cout << endl ;
    fR.push_back(UNSET) ;
  }
  regn_vlist.push_back(nregions) ;
}
void mesh2D_reader_FaceFormat::read_faces( ifstream & inpf ) {
  int nfaces, id, node0, node1, regn0, regn1, dum0 ;
  inpf >> nfaces >> dum0 ;
  // ----
  //cout << "nfaces = " << nfaces << endl ;
  for ( int i=0; i<nfaces; ++i ) {
    inpf >> id >> node0 >> node1 >> regn0 >> regn1 ;
    regn1 = regn1<0 ? UNSET+offset : regn1 ;
    face_vlist.push_back(node0-offset) ;
    face_vlist.push_back(node1-offset) ;
    face_rlist.push_back(regn0-offset) ;
    face_rlist.push_back(regn1-offset) ;
    // ----
    fF.push_back(UNSET) ;
    //cout << id << "  " << node0 << "  " << node1 << "  " << regn0 << "  " << regn1 << endl ;
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
  //cout << "offset = " << offset << endl ;
  //cout << endl ;
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

/////////////////////////////////////////////////////////////

class mesh2D_reader_DurhamFormat : public mesh2D_reader {
  
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
  int        offset    ;  // ==1 fortran & matlab offset
  
  // instantiate work arrays
  vector<int>    face_vlist, face_rlist, regn_vlist, regn_flist, fV, fF, fR ; 
  vector<double> xV, yV ; 

  // methods
  void read_MESH  ( ifstream & inpf ) ;
  void read_OFFSET( ifstream & inpf ) ;
  void read_POINTS( ifstream & inpf ) ;
  void read_EDGES ( ifstream & inpf ) ;
  void read_CELLS_POINTS( ifstream & inpf ) ;
  void read_CELLS_EDGES ( ifstream & inpf ) ;

  void read_file  () ;
  bool read_line  ( ifstream & inpf ) ;

  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV ) { assert(false) ; }
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) { assert(false) ; }
  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR ) {
    read_file() ;
  }

public:
  mesh2D_reader_DurhamFormat( mesh_2Dv & _inp_mesh, string _file_name, int _offset=1 ) : 
    inp_mesh(_inp_mesh), file_name(_file_name), offset(_offset) {}
  ~mesh2D_reader_DurhamFormat() {}
  
  void read_and_build( int mesh_flag ) ;
} ;

void mesh2D_reader_DurhamFormat :: read_and_build( int mesh_flag ) {
  DBGF(begin  -->> mesh2D_reader_DurhamFormat) ;
  
  // input mesh from data
  read_file() ;

  if (false) { // debugging, TBR
    for ( int i=0; i<xV.size(); ++i ) {
      printf("i=%i  xV[%i]=%14.7e  yV[%i]=%14.7e  fV[%i]=%i\n",i,i,xV[i],i,yV[i],i,fV[i]) ;
    }
    for ( int i=0; i<face_vlist.size(); ++i ) {
      printf("i=%i  face_vlist[%i]=%2i  face_rlist[%i]=%2i\n",i,i,face_vlist[i],i,face_rlist[i]) ;
    }
  }

  // start builder
  mesh2Dv_builder mesh_builder( inp_mesh ) ;

  switch ( mesh_flag ) {
  case 0:
    mesh_builder . build_the_mesh( xV,yV, fV, regn_vlist, fR ) ;
    break ;
  case 1:
    mesh_builder . build_the_mesh( xV,yV, fV, face_vlist, face_rlist, fF ) ;
    break ;
  case 2:
    mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, regn_flist, face_vlist, face_rlist, fF ) ;
    break ;
  default:
    assert(false) ;
  }

  inp_mesh.set_mesh_name(string("Mesh")) ;
  DBGF(end of -->> mesh2D_reader_DurhamFormat) ;
}

void mesh2D_reader_DurhamFormat::read_MESH( ifstream & inpf ) {
  int vrs ;
  inpf >> vrs ;
  //cout << " version = " << vrs << endl ;
}

void mesh2D_reader_DurhamFormat::read_OFFSET( ifstream & inpf ) {
  inpf >> offset ;
  // ----
  //cout << "OFFSET " << offset << endl ;
}

void mesh2D_reader_DurhamFormat::read_POINTS( ifstream & inpf ) {
  int nnodes ;
  double x, y ;
  inpf >> nnodes ;
  cout << "nnodes = " << nnodes << endl ;
  for ( int i=0; i<nnodes; ++i ) {
    inpf >> x >> y ;
    xV.push_back(x) ;
    yV.push_back(y) ;
    fV.push_back(UNSET) ;
    // -----
    //cout << std::scientific 
    // << std::setprecision(14) 
    // << "  " << x
    // << "  " << y
    // << endl ;
  }
}
void mesh2D_reader_DurhamFormat::read_CELLS_POINTS( ifstream & inpf ) {
  int nregions, nnodes, inode ;
  inpf >> nregions ;
  // ----
  //cout << "CELL_POINTS " << nregions << endl ;
  for ( int i=0; i<nregions; ++i ) {
    inpf >> nnodes ;
    regn_vlist.push_back(nnodes) ;
    // -----
    //cout << "  " << nnodes << "\t" ;
    // -----
    for ( int j=0; j<nnodes; ++j ) {
      inpf >> inode ;
      regn_vlist.push_back(inode-offset) ;
      // -----
      //cout << inode << "\t" ;
    }
    //cout << endl ;
    fR.push_back(UNSET) ;
  }
  regn_vlist.push_back(nregions) ;
}
void mesh2D_reader_DurhamFormat::read_CELLS_EDGES( ifstream & inpf ) {
  int nregions, nfaces, iface ;
  inpf >> nregions ;
  // ----
  //cout << "CELL_EDGES " << nregions << endl ;
  for ( int i=0; i<nregions; ++i ) {
    inpf >> nfaces ;
    regn_flist.push_back(nfaces) ;
    // -----
    //cout << "  " << nfaces << "\t" ;
    for ( int j=0; j<nfaces; ++j ) {
      inpf >> iface ;
      regn_flist.push_back(iface-offset) ;
      // -----
      //cout << iface << "\t" ;
    }
    //cout << endl ;
    //fR.push_back(UNSET) ;
  }
  regn_flist.push_back(nregions) ;
}
void mesh2D_reader_DurhamFormat::read_EDGES( ifstream & inpf ) {
  int nfaces, id, node0, node1, regn0, regn1, dum0 ;
  inpf >> nfaces ;
  // ----
  cout << "EDGES " << nfaces << endl ;
  for ( int i=0; i<nfaces; ++i ) {
    inpf >> node0 >> node1 >> regn0 >> regn1 ;
    //regn1 = regn1<0 ? UNSET+offset : regn1 ;
    if ( regn0==offset-1 ) {
      regn0 = UNSET+offset ;
      face_vlist.push_back(node1-offset) ;
      face_vlist.push_back(node0-offset) ;
      face_rlist.push_back(regn1-offset) ;
      face_rlist.push_back(regn0-offset) ;
    } else {
      regn1 = regn1==offset-1 ? UNSET+offset : regn1 ;
      face_vlist.push_back(node0-offset) ;
      face_vlist.push_back(node1-offset) ;
      face_rlist.push_back(regn0-offset) ;
      face_rlist.push_back(regn1-offset) ;
    }
    // ----
    fF.push_back(UNSET) ;
    //cout << "  "
    // << node0 << "\t" << node1 << "\t" 
    // << regn0 << "\t" << regn1 << endl ;
  }
  //cout << endl ;
  face_vlist.push_back(nfaces) ;
  face_rlist.push_back(nfaces) ;
}
void mesh2D_reader_DurhamFormat::read_file() {
  string keywd = "" ;
  ifstream inpf( file_name.c_str() ) ;
  if ( !inpf.good() ) { fatal_error(file_name) ; }
  while ( !inpf.eof() ) { read_line( inpf ) ; }
  inpf.close()  ;
}
bool mesh2D_reader_DurhamFormat :: read_line( ifstream & inpf ) {
  string keywd = "" ;
  inpf >> get_comments >> keywd ;
  //PRT(keywd) ;
  bool retval = true ;
  if      ( keywd == "MESH"         ) { read_MESH  (inpf) ; }
  else if ( keywd == "OFFSET"       ) { read_OFFSET(inpf) ; }
  else if ( keywd == "POINTS"       ) { read_POINTS(inpf) ; }
  else if ( keywd == "CELLS_POINTS" ) { read_CELLS_POINTS(inpf) ; }
  else if ( keywd == "CELLS_EDGES"  ) { read_CELLS_EDGES (inpf) ; }
  else if ( keywd == "EDGES"        ) { read_EDGES (inpf) ; }
  else                                { retval = false    ; }
  if ( retval )                       { inpf >> eatline   ; } 
  return retval ;
}

#endif // end of _MESH_2D_GREADER_HH
