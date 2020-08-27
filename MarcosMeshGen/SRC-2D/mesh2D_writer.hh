#include <ctime>

class mesh2Dv_writer {
private:
  mesh_2Dv & mesh ;
  int offset ;
  
public:
  mesh2Dv_writer( mesh_2Dv & _mesh, int _offset=OFFSET  ) : mesh(_mesh), offset(_offset) {}
  ~mesh2Dv_writer() {}

  void write_mesh            ( string fname=string("mesh2D") ) ; 
  void write_mesh_TRIANGLE   ( string fname=string("mesh2D") ) ; 
  void write_mesh_with_marker( string fname=string("mesh2D") ) ; 
  void write_mesh_PT_format  ( string fname=string("mesh2D") ) ; 
  void write_mesh_Suku_format( string fname=string("mesh2D") ) ; 
  void write_mesh_Durham_format( string fname=string("mesh2D"), string comment=string("") ) ; 
  
  int  get_offset()              { return offset ; }
  void set_offset( int _offset ) { offset=_offset ; }
} ;

void mesh2Dv_writer :: write_mesh( string fname ) {
  MSG("start  mesh2Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;
  
  string fnode = fname + string(".node") ;
  ofstream out_node(fnode.c_str()) ;
  out_node << "# *.node file of 2D-mesh in REGN_FACE format " << endl ;
  out_node << "# " << fnode << endl ;
  out_node << mesh.n_vertex() << "  2  0  0" << endl ;
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    out_node << "  " << iV+offset 
	     << "  " << mesh.coords_V(iV,0) 
	     << "  " << mesh.coords_V(iV,1) 
	     << endl ;
  }
  out_node << "# output from mesh2Dv_writer.hh " << endl ;
  out_node.close() ;

  /*
    header line:
    nR nRV 1   [1=dump one region flag]
    
    if nRV=0  =>  variable format:
    [ iR nRV  iV_[0] ... iV_[nRV]   fR ] for iR=0..nR
  */
  int nR = mesh.n_region() ;
  string fele = fname + string(".ele") ;
  ofstream out_ele(fele.c_str()) ;
  out_ele << "# *.ele file of 2D-mesh in REGN_FACE format " << endl ;
  out_ele << "# " << fele << endl ;
  out_ele << nR << "   0   0" << endl ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    vector<int> R_vlist ;
    mesh.get_regn_vrtx( iR, R_vlist ) ;
    int nRV = R_vlist.size() ;
    out_ele << iR + offset << "  " << nRV << "  " ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      out_ele << "  " << R_vlist[ilV]+offset ;
    }
    out_ele << endl ;
  }
  out_ele << "# output from mesh2Dv_writer.hh " << endl ;
  out_ele.close() ;

  MSG("end of mesh2Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush) ;
}

void mesh2Dv_writer :: write_mesh_TRIANGLE( string fname ) {
  MSG("start  mesh2Dv_writer_TRIANGLE_format::write_mesh"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;
  
  string fnode = fname + string(".node") ;
  ofstream out_node(fnode.c_str()) ;
  out_node << "# *.node file of 2D-mesh in TRIANGLE format " << endl ;
  out_node << "# " << fnode << endl ;
  out_node << mesh.n_vertex() << "  2  0  1" << endl ;
  out_node << "# node, x, y, flag " << endl ;
  out_node << "#   flag=0 ==> internal node " << endl ; 
  out_node << "#   flag=1 ==> boundary node with Dirichlet condition" << endl ;
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    out_node << "  " << iV+offset 
	     << "  " << mesh.coords_V(iV,0) 
	     << "  " << mesh.coords_V(iV,1)
	     << "  " << mesh.is_boundary_vrtx(iV) 
	     << endl ;
  }
  out_node << "# output from mesh2Dv_writer.hh " << endl ;
  out_node.close() ;

  /*
    header line:
    nR nRV 1   [1=dump one region flag]
    
    if nRV=0  =>  variable format:
    [ iR nRV  iV_[0] ... iV_[nRV]   fR ] for iR=0..nR
  */
  int nR = mesh.n_region() ;
  string fele = fname + string(".ele") ;
  ofstream out_ele(fele.c_str()) ;
  out_ele << "# *.ele file of 2D-mesh in TRIANGLE format " << endl ;
  out_ele << "# " << fele << endl ;
  out_ele << nR << "   0   0" << endl ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    vector<int> R_vlist ;
    mesh.get_regn_vrtx( iR, R_vlist ) ;
    int nRV = R_vlist.size() ;
    assert( nRV==3 ) ; // only for triangles
    out_ele << iR + offset << "  " << "  " ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      out_ele << "  " << R_vlist[ilV]+offset ;
    }
    out_ele << endl ;
  }
  out_ele << "# output from mesh2Dv_writer.hh " << endl ;
  out_ele.close() ;

  MSG("end of mesh2Dv_writer_TRIANGLE_format::write_mesh"<<endl<<flush) ;
}

// add the marker at the end
void mesh2Dv_writer :: write_mesh_with_marker( string fname ) {
  MSG("start  mesh2Dv_writer::write_mesh_with_marker_format"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;

  string fnode = fname + string(".node") ;
  ofstream out_node(fnode.c_str()) ;
  out_node << "# *.node file of 2D-mesh in REGN_FACE format " << endl ;
  out_node << "# " << fnode << endl ;
  out_node << mesh.n_vertex() << "  2  0  1" << endl ;
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    out_node << "  " << iV+offset 
	     << "  " << mesh.coords_V(iV,0) 
	     << "  " << mesh.coords_V(iV,1) 
	     << "  " << mesh.get_fV  (iV) 
	     << endl ;
  }
  out_node << "# output from mesh2Dv_writer.hh " << endl ;
  out_node.close() ;

  /*
    header line:
    nR nRV 1   [1=dump one region flag]
    
    if nRV=0  =>  variable format:
    [ iR nRV  iV_[0] ... iV_[nRV]   fR ] for iR=0..nR
  */
  int nR = mesh.n_region() ;
  string fele = fname + string(".ele") ;
  ofstream out_ele(fele.c_str()) ;
  out_ele << "# *.ele file of 2D-mesh in REGN_FACE format " << endl ;
  out_ele << "# " << fele << endl ;
  out_ele << nR << "   0   1" << endl ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    vector<int> R_vlist ;
    mesh.get_regn_vrtx( iR, R_vlist ) ;
    int nRV = R_vlist.size() ;
    out_ele << iR + offset << "  " << nRV << "  " ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      out_ele << "  " << R_vlist[ilV]+offset ;
    }
    out_ele << "   " << mesh.get_fR(iR) ;
    out_ele << endl ;
  }
  out_ele << "# output from mesh2Dv_writer.hh " << endl ;
  out_ele.close() ;

  MSG("end of mesh2Dv_writer::write_mesh_with_marker_format"<<endl<<flush) ;
}

// Paulino-Talischi format
void mesh2Dv_writer :: write_mesh_PT_format( string fname ) {
  MSG("start  mesh2Dv_writer_TRIANGLE_format::write_mesh"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;

  ofstream outf(fname.c_str()) ;
  outf << "%HEADER " << endl ;
  outf << "%NODE LIST " << endl ;
  outf << mesh.n_vertex() << endl ;
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    outf << "  " << iV+offset 
	 << "  " << mesh.coords_V(iV,0) 
	 << "  " << mesh.coords_V(iV,1)
	 << endl ;
  }
  /*
    header line:
    nR nRV 1   [1=dump one region flag]
    
    if nRV=0  =>  variable format:
    [ iR nRV  iV_[0] ... iV_[nRV]   fR ] for iR=0..nR
  */
  int nR = mesh.n_region() ;
  outf << "%ELEMENT CONNECTIVITY" << endl ;
  outf << nR << endl ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    vector<int> R_vlist ;
    mesh.get_regn_vrtx( iR, R_vlist ) ;
    int nRV = R_vlist.size() ;
    outf << iR + offset << "  " << nRV ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      outf << "  " << R_vlist[ilV]+offset ;
    }
    outf << endl ;
  }
  outf << "%END" << endl ;
  outf.close() ;

  MSG("end of mesh2Dv_writer_TRIANGLE_format::write_mesh"<<endl<<flush) ;
}

// Sukumar format
void mesh2Dv_writer :: write_mesh_Suku_format( string fname ) {
  MSG("start  mesh2Dv_writer_Sukumar_format::write_mesh"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;

  // offset is always 1 (for Matlab)
  set_offset( 1 ) ;

  { // write node files 
    string fnode = fname + string(".node") ;
    ofstream outf(fnode.c_str()) ;
    for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
      int bnd_flag = mesh.is_boundary_vrtx(iV) ? 1 : 0 ;
      outf << std::setw(20) 
	   << std::setprecision(14) 
	   << "  " << mesh.coords_V(iV,0) 
	   << "  " << mesh.coords_V(iV,1)
	   << std::setw(4)
	   << "  "  << bnd_flag 
	   << endl ;
    }
    outf.close() ;
  }

  { // write ele files
    string fele = fname + string(".ele") ;
    ofstream outf(fele.c_str()) ;
    int nR = mesh.n_region() ;
    for ( int iR=0 ; iR<nR ; ++iR ) {
      vector<int> R_vlist ;
      mesh.get_regn_vrtx( iR, R_vlist ) ;
      int nRV = R_vlist.size() ;
      outf << nRV ;
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	outf << "  " << R_vlist[ilV]+offset ;
      }
      outf << endl ;
    }
    outf.close() ;
  }

  MSG("end of mesh2Dv_writer_Sukumar_format::write_mesh"<<endl<<flush) ;
}

// Durham format
void mesh2Dv_writer :: write_mesh_Durham_format( string fname, string comment ) {
  MSG("start  mesh2Dv_writer_Durham_format::write_mesh"<<endl<<flush) ;
  MSG("write \""<<fname<<"\""<<endl<<flush) ;

  fname += string(".mesh") ;
  ofstream outf(fname.c_str()) ;
    
  // offset is always 1 (for Matlab)
  set_offset( 1 ) ;

  { // comments
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    outf << "# mesh in Durham's format" << endl ;
    outf << "# output generated on " << dt ;
    outf << "# " << comment << endl ;
    outf << "# " << endl ;
  }

  { // write MESH record
    outf << "MESH 1.00" << endl ;
    outf << "#" << endl ;
  }
  
  { // write OFFSET record
    outf << "OFFSET " << offset << endl ;
    outf << "#" << endl ;
  }
  
  { // write POINTS records
    int nV = mesh.n_vertex() ;
    outf << "POINTS " << nV << endl ;
    for ( int iV=0 ; iV<nV ; ++iV ) {
      outf << std::scientific 
	   << std::setprecision(14) 
	   << "  " << mesh.coords_V(iV,0) 
	   << "  " << mesh.coords_V(iV,1)
	   << endl ;
    }
    outf << "#" << endl ;
  }

  { // write CELLS_POINTS files
    int nR = mesh.n_region() ;
    outf << "CELLS_POINTS " << nR << endl ;
    for ( int iR=0 ; iR<nR ; ++iR ) {
      vector<int> regn_vlist ;
      mesh.get_regn_vrtx( iR, regn_vlist ) ;
      int nRV = regn_vlist.size() ;
      outf << "  " << nRV << "\t" ;
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	outf << regn_vlist[ilV]+offset << "\t" ;
      }
      outf << endl ;
    }
    outf << "#" << endl ;
  }

  { // write CELLS_EDGES files
    int nR = mesh.n_region() ;
    outf << "CELLS_EDGES " << nR << endl ;
    for ( int iR=0 ; iR<nR ; ++iR ) {
      vector<int> regn_flist ;
      mesh.get_regn_face( iR, regn_flist ) ;
      int nRF = regn_flist.size() ;
      outf << "  " << nRF << "\t" ;
      for ( int ilF=0; ilF<nRF; ++ilF ) {
	outf << regn_flist[ilF]+offset << "\t" ;
      }
      outf << endl ;
    }
    outf << "#" << endl ;
  }

  { // write EDGES files
    int nF = mesh.n_face() ;
    outf << "EDGES " << nF << endl ;
    for ( int iF=0 ; iF<nF ; ++iF ) {
      int iR1 = mesh.is_boundary_face(iF) ? offset-1 : mesh.face_regn(iF,1)+offset ;
      outf << "  "
	   << mesh.face_vrtx(iF,0)+offset << "\t"
	   << mesh.face_vrtx(iF,1)+offset << "\t"
	   << mesh.face_regn(iF,0)+offset << "\t" 
	   << iR1 << endl ;
    }
    outf << "#" << endl ;
  }

  // final, close the output file
  outf.close() ;

  MSG("end of mesh2Dv_writer_Durham_format::write_mesh"<<endl<<flush) ;
}
