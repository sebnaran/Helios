void mesh2D_reader_DurhamFormat :: read_the_mesh( vector<double> & xV, 
						  vector<double> & yV,  
						  vector<int>    & fV, 
						  vector<int>    & regn_vlist, 
						  vector<int>    & fR ) {
  // read DATA file
  string data_file_name = file_name + ".data" ;
  ifstream data_file( data_file_name.c_str() ) ;
  if ( data_file.good() ) {
    read_file( xV, yV, fV, regn_vlist, fR ) ;
  } else {
    error_message(data_file_name) ;
  }
  data_file.close() ;
}
void mesh2D_reader_DurhamFormat :: read_file( string _fname="" ) {
  string inp_fname = _fname=="" ? fname : _fname ;
  ifstream inpf( inp_fname.c_str() ) ;
  if ( !inpf.good() ) { fatal_error(inp_fname) ; }
    while ( !inpf.eof() ) { read_line( inpf ) ; }
    inpf.close()  ;
}

bool mesh2D_reader_DurhamFormat :: read_line( string keywd, ifstream & inpf ) {
  bool retval = true ;
  if      ( keywd == "MESH"   ) { read_MESH  (inpf) ; }
  else if ( keywd == "POINTS" ) { read_POINTS(inpf) ; }
  else if ( keywd == "CELLS"  ) { read_CELLS (inpf) ; }
  else if ( keywd == "EDGES"  ) { read_EDGES (inpf) ; }
  else                          { retval = false    ; }
  if ( retval )                 { inpf >> eatline   ; } 
  return retval ;
} 

void mesh2D_reader_DurhamFormat :: read_vrtx_data ( ifstream & input_file, 
						    vector<double> & xV, 
						    vector<double> & yV,
						    vector<int>    & fV ) {
  MSG("begin mesh2D_reader_DurhamFormat::read_vrtx_data"<<endl) ;
  
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
  MSG("end of mesh_reader_DurhamFormat::read_vrtx_data"<<endl) ;
}
void mesh2D_reader_DurhamFormat :: read_regn_data( ifstream & input_file, 
						     vector<int> & regn_vlist, 
						     vector<int> & fR ) {
  MSG("begin mesh2D_reader_DurhamFormat::read_regn_data"<<endl) ;
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
  MSG("end of mesh2D_reader_DurhamFormat::read_regn_data"<<endl) ;
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
