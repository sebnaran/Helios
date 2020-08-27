#ifndef _MESH_2D_PSOUT_HH
#define _MESH_2D_PSOUT_HH

/**

VERY SIMPLE SEQUENCE FOR DRAWING MESHES

// define the string containing mesh name from mesh "current_mesh"
string mesh_name = current_mesh.get_mesh_name() ;

// define a ps driver for the file "mesh_name"
PS_Driver ps_driver( mesh_name ) ;

// define a psout object for drawing items from mesh "current_mesh"
mesh2D_psout<common_2D> mesh2D_psout( current_mesh, ps_driver ) ;

// no numbering of mesh items (default is true)
bool draw_numbers(false) ;

// call a drawing method
mesh2D_psout.draw_mesh_0( mesh_name, draw_numbers ) ;

**/

class PS_Driver {
  
private:
  // private pars
  double ps_xo    ;
  double ps_yo    ;
  double ps_width ;
  
  // mesh bounding box pars
  double ps_xmin  ;
  double ps_ymin  ;
  double ps_delta ;
  
  // output stream file
  ofstream psout ;
  bool status ;
  
  // interface for postscript position vector (X,Y)
  double X( double x ) const { return ps_xo + ps_width * ( x-ps_xmin )/ps_delta ; }
  double Y( double y ) const { return ps_yo + ps_width * ( y-ps_ymin )/ps_delta ; }
  
public: 
  // CTOR/DTOR
  PS_Driver() : status(false) {}
  PS_Driver( string fname ) : status(false) {
    setup(fname) ;
  }
  ~PS_Driver() { close() ; }
  
  // display methods  
  void setup( string fname ) {
    // if already opened, close it
    if ( status ) { close() ; }
    else          { status=bool(true) ; }
    // open file
    string inp_file( fname+string(".ps") ) ;

    psout.open( inp_file.c_str() ) ;
    if ( psout.good() ) { 
      status = bool(true) ;
      open_postscript_file() ; 
    }
    else { 
      cerr << "Cannot open for write file: ``" << inp_file << "''" << endl ; 
    }
    //setup pars
    ps_xo    =  50 ;
    ps_yo    = 200 ;
    ps_width = 500 ;
  }

  void open_postscript_file() ;
  void close() { 
    status = bool(false) ;
    clear() ; 
    psout.close() ; 
  }
  void clear() { psout << endl << "showpage" << endl ; }
  void setup_bbox( double xmin, double ymin, double xmax, double ymax ) {
    ps_xmin  = xmin ;
    ps_ymin  = ymin ;
    ps_delta = max(xmax-xmin, ymax-ymin) ;
  }
  void draw_ellp( double x0, double y0, double width=-1. ) {
    width = width>0. ? width : ps_width/100. ;
    psout << width << " " << width << " "
	  << X(x0) << " " << Y(y0) << " ELLF" << endl ;
  }
  void draw_line( double xstr, double ystr, double xend, double yend ) {
    psout << X(xstr) << " " << Y(ystr) << " "
	  << X(xend) << " " << Y(yend) << " DL" << endl ; 
  }
  void draw_rect( double x0, double y0, double x1, double y1 ) {
    psout << X(x0) << " " << Y(y0) << " "
	  << X(x1) << " " << Y(y1) << " DR " << endl ; 
  }
  void fill_rect( double x0, double y0, double x1, double y1 ) {
    psout << X(x0) << " " << Y(y0) << " "
	  << X(x1) << " " << Y(y1) << " FR " << endl ; 
  }
  void draw_tria( double xa, double ya, double xb, double yb, double xc, double yc ) {
    psout << X(xa) << "  " << Y(ya) << "  "
	  << X(xb) << "  " << Y(yb) << "  "
	  << X(xc) << "  " << Y(yc) << "  DT " << endl ; 
  }
  void fill_tria( double xa, double ya, double xb, double yb, double xc, double yc ) {
    psout << X(xa) << "  " << Y(ya) << "  "
	  << X(xb) << "  " << Y(yb) << "  "
	  << X(xc) << "  " << Y(yc) << "  FT " << endl ; 
  }
  void dump_string( string str ) { psout << str << endl ; }
  void draw_string( string str, double xs, double ys ) {
    psout << " (" << str << ") " << X(xs) << " " << Y(ys) << " VLS" << endl ;
  }
  void draw_string( string str, int n, double xs, double ys ) {
    psout << " (" << str << " " << n << ") " << X(xs) << " " << Y(ys) << " VLS" << endl ;
  }
  void draw_number( int n, double xn, double yn ) {
    psout << " (" << n << ") " << X(xn) << "  " << Y(yn) << "  VLS " << endl ;
  }
  // setup colors
  void set_rgb_colour( string colour ) { psout << " " << colour << " " << endl ; }
  void set_line_type ( string ltype  ) { psout << " " << ltype  << " " << endl ; }
  void set_font      ( string font )   { 
    psout << string(" set") + font + string("LblFont") << endl ; 
  }
} ;

void PS_Driver::open_postscript_file() {
  psout << "%!PS"                                            << endl
        << "%%BoundingBox:    0   125   600   850"           << endl
        << endl
        << "%--------------------------------------------"   << endl
        << "%----------- PS Cmd Definitions--------------"   << endl
        << "%--------------------------------------------"   << endl
        << endl
        << "/DL      { moveto lineto stroke } def"           << endl
        << "/TT      { translate } bind def"                 << endl
        << "/SS      { scale     } bind def"                 << endl
        << "/RESTORE { cmtx setmatrix } def"                 << endl
        << "/SAVE    { /cmtx matrix currentmatrix def } def" << endl
        << "/ELL     { gsave SAVE TT SS newpath"             << endl
        << "           0 0 1 0 360 arc RESTORE closepath"    << endl
        << "           grestore } def"                       << endl
        << "/ELLF    % drawellipse -- stack: sx,sy,x0,y0"    << endl
        << "         { gsave SAVE TT SS newpath"             << endl
        << "           0 0 1 0 360 arc RESTORE closepath"    << endl
        << "           fill grestore } def"                  << endl
        << "/DR       % drawrect -- stack : x0,y0,xl,yl"     << endl 
        << "          { /ywidth exch def"                    << endl 
        << "            /xwidth exch def"                    << endl      
        << "            newpath moveto"                      << endl
        << "            0 ywidth rlineto"                    << endl
        << "            xwidth 0 rlineto"                    << endl
        << "            0 ywidth neg rlineto"                << endl
        << "            closepath"                           << endl
        << "            } def"                               << endl
        << "/FR       % fillrect -- stack : x0,y0,xl,yl"     << endl
        << "          { DR fill stroke } def"                << endl 
        << "/DT       % stack : x0,y0,x1,y1,x2,y2"           << endl 
        << "          { /p3y exch def"                       << endl 
        << "            /p3x exch def"                       << endl 
        << "            /p2y exch def"                       << endl 
        << "            /p2x exch def"                       << endl 
        << "            newpath moveto"                      << endl 
        << "            p2x p2y lineto"                      << endl 
        << "            p3x p3y lineto"                      << endl 
        << "            closepath fill stroke"               << endl 
        << "            } def"                               << endl 
        << "/FT       % stack : x0,y0,x1,y1,x2,y2"           << endl 
        << "          { DT fill } def"                 << endl 
        << endl
        << "%--------------------------------------------"   << endl
        << "%--------------Font Definitions--------------"   << endl
        << "%--------------------------------------------"   << endl
        << endl
        << "/findheight   % Find the font height"            << endl
        << "{ gsave"                                         << endl
        << "  newpath"                                       << endl
        << "  0 0 moveto"                                    << endl
        << "  (X) true charpath"                             << endl
        << "  flattenpath"                                   << endl
        << "  pathbbox /capheight exch def pop pop pop"      << endl
        << "  grestore"                                      << endl
        << "} def"                                           << endl
        << endl
        << "% Load ISOLatin font for Times-Roman"            << endl
        << "/Times-Roman findfont"                           << endl
        << "dup length dict begin"                           << endl
        << "  {1 index /FID ne {def} {pop pop} ifelse} forall" << endl
        << "  /Encoding ISOLatin1Encoding def"               << endl
        << "  currentdict"                                   << endl
        << "end"                                             << endl
        << "/Times-Roman-ISOLatin1 exch definefont pop"      << endl
        << endl
        << "% Load ISOLatin font for Times-Bold"             << endl
        << "/Times-Bold findfont"                            << endl
        << "dup length dict begin"                           << endl
        << "  {1 index /FID ne {def} {pop pop} ifelse} forall" << endl
        << "  /Encoding ISOLatin1Encoding def"               << endl
        << "  currentdict"                                   << endl
        << "end"                                             << endl
        << "/Times-Bold-ISOLatin1 exch definefont pop"       << endl
        << endl
        << "%------------------------------------------"     << endl
        << "%- Vrtx/Edge/Poly Label Macro Definitions -"     << endl
        << "%------------------------------------------"     << endl
        << endl
        << "% Define font sizes"                             << endl
        << "/VrtxLblFontSize 12.00 def"                      << endl
        << "/EdgeLblFontSize 12.00 def"                      << endl
        << "/PolyLblFontSize 12.00 def"                      << endl
        << "/MeshLblFontSize 24.00 def"                      << endl
        << endl
        << "% Define fonts"                                  << endl
        << "/VrtxLblFont /Times-Bold-ISOLatin1 findfont"     << endl
        << "VrtxLblFontSize scalefont def"                   << endl
        << "/EdgeLblFont /Times-Bold-ISOLatin1 findfont"     << endl
        << "EdgeLblFontSize scalefont def"                   << endl
        << "/PolyLblFont /Times-Bold-ISOLatin1 findfont"     << endl
        << "PolyLblFontSize scalefont def"                   << endl
        << "/MeshLblFont /Times-Bold-ISOLatin1 findfont"     << endl
        << "MeshLblFontSize scalefont def"                   << endl
        << endl
        << "% Define font locations/offsets"                 << endl
        << "/VrtxLblOffset {VrtxLblFontSize -1 div"          << endl
        << "                VrtxLblFontSize -1 div rmoveto } def" << endl
        << "/EdgeLblOffset {EdgeLblFontSize -2 div"          << endl
        << "                EdgeLblFontSize -3 div rmoveto } def" << endl
        << "/PolyLblOffset {PolyLblFontSize -3 div"          << endl
        << "                PolyLblFontSize    neg rmoveto } def" << endl
        << "/MeshLblOffset {MeshLblFontSize -3 div"          << endl
        << "                MeshLblFontSize    neg rmoveto } def" << endl
        << endl
        << "% Define font definition macros" << endl
        << "/setVrtxLblFont {VrtxLblFont setfont findheight} def" << endl
        << "/setEdgeLblFont {EdgeLblFont setfont findheight} def" << endl
        << "/setPolyLblFont {PolyLblFont setfont findheight} def" << endl
        << "/setMeshLblFont {MeshLblFont setfont findheight} def" << endl
        << endl 
        << "% Define Vrtx/Edge/Poly show commands" << endl 
        << "/VLS { moveto VrtxLblOffset show } def" << endl 
        << "/ELS { moveto EdgeLblOffset show } def" << endl 
        << "/PLS { moveto PolyLblOffset show } def" << endl 
        << "/TLS { moveto MeshLblOffset show } def" << endl 
        << endl 
        << "%------------------------------------------" << endl
        << "%----------Color Table Definitions---------" << endl
        << "%------------------------------------------" << endl
        << endl
        << "% set colours" << endl
        << "/BLACK     { 0.00 0.00 0.00 setrgbcolor } def" << endl
        << "/RED       { 1.00 0.00 0.00 setrgbcolor } def" << endl
        << "/GREEN     { 0.00 1.00 0.00 setrgbcolor } def" << endl
        << "/MAGENTA   { 1.00 0.00 1.00 setrgbcolor } def" << endl
        << "/BROWN     { 0.60 0.35 0.25 setrgbcolor } def" << endl
        << "/BLUE      { 0.00 0.00 1.00 setrgbcolor } def" << endl
	<< "/LIGHTFILL { 0.90 0.90 0.90 setrgbcolor } def" << endl
        << endl 
        << "%------------------------------------------" << endl
        << "%----------LType Table Definitions---------" << endl
        << "%------------------------------------------" << endl
        << endl
        << "/SOLID {   [] 0 setdash } def" << endl
        << "/DASH1 {  [5] 0 setdash } def" << endl
        << "/DASH2 { [10] 0 setdash } def" << endl
        << endl
        << "/MLT3  {[1 2] 1 setdash } def   % Dotted Line" << endl
        << "/GRAY00 { 1.00 setgray } def" << endl
        << "/GRAY01 { 1.05 setgray } def" << endl
        << "/GRAY02 { 1.10 setgray } def" << endl
        << "/GRAY03 { 1.15 setgray } def" << endl
        << "/GRAY04 { 1.20 setgray } def" << endl
        << "/GRAY05 { 1.25 setgray } def" << endl
        << "/GRAY06 { 1.30 setgray } def" << endl
        << "/GRAY07 { 1.35 setgray } def" << endl
        << "/GRAY08 { 1.40 setgray } def" << endl
        << "/GRAY09 { 1.45 setgray } def" << endl
        << "/GRAY10 { 1.50 setgray } def" << endl
        << "/GRAY11 { 1.55 setgray } def" << endl
        << "/GRAY12 { 1.60 setgray } def" << endl
        << "/GRAY13 { 1.65 setgray } def" << endl
        << "/GRAY14 { 1.70 setgray } def" << endl
        << "/GRAY15 { 1.75 setgray } def" << endl
        << "/GRAY16 { 1.80 setgray } def" << endl
        << "/GRAY17 { 1.85 setgray } def" << endl
        << "/GRAY18 { 1.90 setgray } def" << endl
        << "/GRAY19 { 1.95 setgray } def" << endl 
	<< endl
        << "newpath"                      << endl 
	<< endl ;
}

#if 0
// graphic driver for mesh classes derived from p2_mesh 
class mesh2D_psout {
  
private:
  mesh_2Dv  & mesh ;
  PS_Driver & psd  ;

  int offset ;

public:
  mesh2D_psout( mesh_2Dv & _mesh, PS_Driver & _psd ) : 
    mesh(_mesh), psd(_psd), offset(0) {}
  ~mesh2D_psout() {}

  void setup_bb_mesh() const ;
  
  void draw_ET_graph() const ;
  void draw_EE_graph() const ;
  void draw_TT_graph() const ;

  void draw_bbox    () const ;
  void draw_boundary() const ;

  void dump_string  ( string str ) const ;

  // draw Vrtx/Edge/Regn position numbers
  void draw_vrtx_pos() const ;
  void draw_edge_pos  () const ;
  void draw_regn_pos  () const ;

  void draw_primal_mesh() const ;
  void draw_median_dual_mesh() const ;
  void draw_center_dual_mesh() const ;

  void draw_mesh_0( string mesh_label_str, bool draw_numbers=true ) const ;
  void draw_mesh_1( string mesh_label_str, bool draw_numbers=true ) const ;
  void draw_mesh  ( string mesh_label_str, bool draw_numbers=true ) const ;

  // offset number; can be changed by mesh (e.g. to start fortran-like numbering from 1)
  void set_offset( int _offset ) { offset = _offset ; } 
  int const get_offset() const { return offset ; }

  void draw_vrtx_bullet    ( int size=0, string col_str=string("BLACK") ) const ;
  void draw_face_bullet    ( int size=0, string col_str=string("RED")   ) const ;
  void draw_regn_bullet    ( int size=0, string col_str=string("BLUE")  ) const ;

  void draw_bnd_vrtx_bullet( int size=0, string col_str=string("BLACK") ) const ;
  void draw_int_vrtx_bullet( int size=0, string col_str=string("BLACK") ) const ;

  void draw_bnd_face_bullet( int size=0, string col_str=string("RED")   ) const ;
  void draw_int_face_bullet( int size=0, string col_str=string("RED")   ) const ;
} ;

void mesh2D_psout::dump_string( string str ) const { 
  psd.dump_string( str ) ; 
}
void mesh2D_psout::setup_bb_mesh() const {
  double xmin, xmax, ymin, ymax ;
  mesh.bbox(xmin, ymin, xmax, ymax) ;
  psd.setup_bbox( xmin, ymin, xmax, ymax ) ;
}
void mesh2D_psout::draw_primal_mesh() const {
  DBGF(BBB-00) ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    int    iV0   = mesh.face_vrtx(iF,0) ;
    int    iV1   = mesh.face_vrtx(iF,1) ;
    double xV[2] = { mesh.coords_V(iV0,0), mesh.coords_V(iV1,0) } ;
    double yV[2] = { mesh.coords_V(iV0,1), mesh.coords_V(iV1,1) } ;
    psd.draw_line( xV[0], yV[0], xV[1], yV[1] ) ;
  }
}
void mesh2D_psout::draw_median_dual_mesh() const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double x0 = mesh.coords_R(iR,0) ;
    double y0 = mesh.coords_R(iR,1) ;
    vector<int> flist ;
    mesh.get_regn_face( iR, flist ) ;
    int nRF = flist.size() ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = flist[ilF] ;
      double x1 = mesh.coords_F(iF,0) ;
      double y1 = mesh.coords_F(iF,1) ;
      psd.draw_line( x0, y0, x1, y1 ) ;
    }
  }
}
void mesh2D_psout::draw_center_dual_mesh() const {
  int nF = mesh.n_face () ;
  for ( int iF=0; iF<nF; ++iF ) {
    
    int    iR0 = mesh.face_regn(iF,0) ;
    double xR0 = mesh.coords_R(iR0,0) ;
    double yR0 = mesh.coords_R(iR0,1) ;
    
    if ( mesh.is_boundary_face(iF) ) {
      double xbF = mesh.coords_F(iF,0) ; ;
      double ybF = mesh.coords_F(iF,1) ; ;
      psd.draw_line( xR0, yR0, xbF, ybF ) ;
    } else {
      int    iR1 = mesh.face_regn(iF,1) ;
      double xR1 = mesh.coords_R(iR1,0) ;
      double yR1 = mesh.coords_R(iR1,1) ;
      psd.draw_line( xR0, yR0, xR1, yR1 ) ;
    }
  }
}
void mesh2D_psout::draw_EE_graph () const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {

    vector<int> flist ;
    mesh.get_regn_face(iR,flist) ;

    int nRF = flist.size() ;
    for ( int il0=0; il0<nRF; ++il0 ) {
      int    il1   = (il0+1)%nRF ;
      int    iF0   = flist[il0] ;
      int    iF1   = flist[il1] ;
      double xF[2] = { mesh.coords_F(iF0,0), mesh.coords_F(iF1,0) } ; 
      double yF[2] = { mesh.coords_F(iF0,1), mesh.coords_F(iF1,1) } ; 
      psd.draw_line( xF[0], yF[0], xF[1], yF[1] ) ;
    }
  }
}
void mesh2D_psout::draw_ET_graph() const {
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    double xFm = mesh.coords_F(iF,0) ;
    double yFm = mesh.coords_F(iF,1) ;
    int    iR0 = mesh.face_regn(iF,0) ;
    double xR0 = mesh.coords_R(iR0,0) ;
    double yR0 = mesh.coords_R(iR0,1) ;
    psd.draw_line( xFm, yFm, xR0, yR0 ) ;
    if ( mesh.is_internal_face(iF) ) {
      int    iR1 = mesh.face_regn(iF,1) ;
      double xR1 = mesh.coords_R(iR1,0) ;
      double yR1 = mesh.coords_R(iR1,1) ;
      psd.draw_line( xFm, yFm, xR1, yR1 ) ;
    }
  }
}
void mesh2D_psout::draw_TT_graph() const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    
    double x0 = mesh.coords_R(iR,0) ;
    double y0 = mesh.coords_R(iR,1) ;

    vector<int> flist ;
    mesh.get_regn_face(iR,flist) ;

    int nRF = flist.size() ;
    for ( int il=0; il<nRF; ++il ) {
      int    iF = flist[il] ;
      double x1 = mesh.coords_F(iF,0) ;
      double y1 = mesh.coords_F(iF,1) ;
      psd.draw_line( x0, y0, x1, y1 ) ;
    }
  }
}
void mesh2D_psout::draw_bbox() const {
  double x0, y0, x1, y1 ;
  mesh.bbox( x0, y0, x1, y1 ) ;
  //psd.draw_rect( x0, y0, x1, y1 ) ;
  psd.draw_line( x0, y0, x1, y0 ) ;
  psd.draw_line( x1, y0, x1, y1 ) ;
  psd.draw_line( x1, y1, x0, y1 ) ;
  psd.draw_line( x0, y1, x0, y0 ) ;
}
void mesh2D_psout::draw_boundary() const {
  for ( int ilF=0; ilF<mesh.n_bface(); ++ilF ) {
    int ibF = mesh.get_bnd_face(ilF) ;
    int iV0 = mesh.face_vrtx(ibF,0) ;
    int iV1 = mesh.face_vrtx(ibF,1) ;
    double xV[2] = { mesh.coords_V(iV0,0), mesh.coords_V(iV1,0) } ;
    double yV[2] = { mesh.coords_V(iV0,1), mesh.coords_V(iV1,1) } ;
    psd.draw_line( xV[0], yV[0], xV[1], yV[1] ) ;
  }
}
//draw_vrtx_pos
void mesh2D_psout::draw_vrtx_pos() const  {
  psd.dump_string( "% ----- vrtx numbers" ) ;
  psd.set_rgb_colour( "BLACK" ) ;
  psd.set_font( "Vrtx" ) ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xV = mesh.coords_V(iV,0) ;
    double yV = mesh.coords_V(iV,1) ;
    psd.draw_number( iV+get_offset(), xV, yV ) ;
  }
}
//draw_edge_pos
void mesh2D_psout::draw_edge_pos() const {
  psd.dump_string( "% ----- edge numbers" ) ;
  psd.set_rgb_colour( "RED" ) ;
  psd.set_font( "Edge" ) ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    double xF = mesh.coords_F(iF,0) ;
    double yF = mesh.coords_F(iF,1) ;
    psd.draw_number( iF+get_offset(), xF, yF ) ;
  }
}
//draw_regn_pos
void mesh2D_psout::draw_regn_pos() const {
  psd.dump_string( "% ----- regn numbers" ) ;
  psd.set_rgb_colour( "BLUE" ) ;
  psd.set_font( "Poly" ) ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double xR = mesh.coords_R(iR,0) ;
    double yR = mesh.coords_R(iR,1) ;
    psd.draw_number( iR+get_offset(), xR, yR ) ;
  }
}
//draw_regn_bullet
void mesh2D_psout::draw_regn_bullet( int size, string col_str ) const {
  psd.dump_string( "% ----- regn bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<=0 ) { size=1 ; }
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double xR = mesh.coords_R(iR,0) ;
    double yR = mesh.coords_R(iR,1) ;
    // MODIFIED HMM position -------------------------
    if ( false ) {
      vector<int> R_vlist ;
      mesh.get_regn_vrtx( iR, R_vlist ) ;
      int nRV = R_vlist.size() ;
      double xRc = 0. ; 
      double yRc = 0. ;
      double sum = 0.; 
      for ( int ilV=0; ilV<3; ++ilV ) {
	int iV = R_vlist[ilV] ;
	xRc += mesh.coords_V( iV, 0 ) ;
	yRc += mesh.coords_V( iV, 1 ) ;
      sum += 1. ;
      }
      xRc /= sum ;
      yRc /= sum ;
      psd.draw_ellp( xRc, yRc, size ) ;
    } else {
      psd.draw_ellp( xR, yR, size ) ;
    }
  }
}
//draw_vrtx_pos 
void mesh2D_psout::draw_int_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<=0 ) { size=1 ; }
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    if ( mesh.is_internal_vrtx(iV) ) {
      double xV = mesh.coords_V(iV,0) ;
      double yV = mesh.coords_V(iV,1) ;
      psd.draw_ellp( xV, yV, size ) ;
    }
  }
}   
void mesh2D_psout::draw_bnd_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<=0 ) { size=1 ; }
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    if ( mesh.is_boundary_vrtx(iV) ) {
      double xV = mesh.coords_V(iV,0) ;
      double yV = mesh.coords_V(iV,1) ;
      psd.draw_ellp( xV, yV, size ) ;
    }
  }
} 
void mesh2D_psout::draw_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<=0 ) { size=1 ; }
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xV = mesh.coords_V(iV,0) ;
    double yV = mesh.coords_V(iV,1) ;
    psd.draw_ellp( xV, yV, size ) ;
  }
}
//draw_vrtx_pos
void mesh2D_psout::draw_int_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<0 ) { size=1 ; }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    if ( mesh.is_internal_face(iF) ) {
      double xF = mesh.coords_F(iF,0) ;
      double yF = mesh.coords_F(iF,1) ;
      psd.draw_ellp( xF, yF, size ) ;
    }
  }
}
void mesh2D_psout::draw_bnd_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<0 ) { size=1 ; }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    if ( mesh.is_boundary_face(iF) ) {
      double xF = mesh.coords_F(iF,0) ;
      double yF = mesh.coords_F(iF,1) ;
      psd.draw_ellp( xF, yF, size ) ;
    }
  }
}
void mesh2D_psout::draw_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size<0 ) { size=1 ; }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    double xF = mesh.coords_F(iF,0) ;
    double yF = mesh.coords_F(iF,1) ;
    psd.draw_ellp( xF, yF, size ) ;
  }
}
void mesh2D_psout::draw_mesh_0( string mesh_label_str, bool draw_numbers ) const {
  setup_bb_mesh() ;                  // init bounding box

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ; 
  draw_primal_mesh() ;

  if ( draw_numbers ) {
    draw_vrtx_pos () ;              // draw vrtx/edge/regn numbers
    draw_edge_pos   () ;
    draw_regn_pos   () ;
  }

  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;
}
void mesh2D_psout::draw_mesh_1( string mesh_label_str, bool draw_numbers ) const {
  setup_bb_mesh() ;                 // init bounding box
  draw_bbox() ;                     // draw bounding box

  psd.set_rgb_colour( "BROWN" ) ;   // draw dual mesh
  psd.set_line_type ( "DASH1" ) ; 
  draw_median_dual_mesh() ;

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ;
  draw_primal_mesh() ;
  
  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
    draw_edge_pos() ;
    draw_regn_pos() ;
  }

  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;
}
void mesh2D_psout::draw_mesh( string mesh_label_str, bool draw_numbers ) const {
  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;

  setup_bb_mesh() ;                 // init bounding box
  //draw_bbox() ;                   // draw bounding box

#if 0
  psd.set_rgb_colour( "BLACK" ) ;   // draw dual mesh
  psd.set_line_type ( "DASH1" ) ;
  draw_center_dual_mesh() ;

  draw_regn_bullet    (3) ;
  draw_bnd_vrtx_bullet(3) ;
  draw_bnd_face_bullet(3) ;
#endif

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ; 
  draw_primal_mesh() ;

  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
    draw_edge_pos() ;
    draw_regn_pos() ;
 }
}
#endif

// graphics driver for mesh classes derived from p2_mesh 
class mesh2D_psout {
  
private:
  mesh_2Dv  & mesh ;
  PS_Driver & psd  ;

  int offset ;

public:
  mesh2D_psout( mesh_2Dv & _mesh, PS_Driver & _psd ) : 
    mesh(_mesh), psd(_psd), offset(0) {}
  ~mesh2D_psout() {}

  void setup_bb_mesh() const ;
  
  void draw_ET_graph() const ;
  void draw_EE_graph() const ;
  void draw_TT_graph() const ;

  void draw_bbox    () const ;
  void draw_boundary() const ;

  void draw_filled_regions( vector<bool> & regn_flag_list ) const ;

  void dump_string  ( string str ) const ;

  // draw Vrtx/Edge/Regn position numbers
  void draw_vrtx_pos() const ;
  void draw_edge_pos  () const ;
  void draw_regn_pos  () const ;

  void draw_primal_mesh() const ;
  void draw_median_dual_mesh() const ;
  void draw_center_dual_mesh() const ;

  void draw_mesh_0( string mesh_label_str, bool draw_numbers=true ) const ;
  void draw_mesh_1( string mesh_label_str, bool draw_numbers=true ) const ;
  void draw_mesh  ( string mesh_label_str, bool draw_numbers=true ) const ;
  void draw_mesh  ( string mesh_label_str, vector<bool> regn_flag_list,  
		    bool draw_numbers=true ) const ;

  void draw_meshless( string mesh_label_str, bool draw_numbers=true, bool draw_meshless=true ) const ;

  // offset number; can be changed by mesh (e.g. to start fortran-like numbering from 1)
  void set_offset( int _offset ) { offset = _offset ; } 
  int const get_offset() const { return offset ; }

  void draw_vrtx_bullet    ( int size=0, string col_str=string("BLACK") ) const ;
  void draw_face_bullet    ( int size=0, string col_str=string("RED")   ) const ;
  void draw_regn_bullet    ( int size=0, string col_str=string("BLUE")  ) const ;

  void draw_bnd_vrtx_bullet( int size=0, string col_str=string("BLACK") ) const ;
  void draw_int_vrtx_bullet( int size=0, string col_str=string("BLACK") ) const ;

  void draw_bnd_face_bullet( int size=0, string col_str=string("RED")   ) const ;
  void draw_int_face_bullet( int size=0, string col_str=string("RED")   ) const ;
} ;

void mesh2D_psout::dump_string( string str ) const { 
  psd.dump_string( str ) ; 
}
void mesh2D_psout::setup_bb_mesh() const {
  double xmin, xmax, ymin, ymax ;
  mesh.bbox(xmin, ymin, xmax, ymax) ;
  psd.setup_bbox( xmin, ymin, xmax, ymax ) ;
}
void mesh2D_psout::draw_primal_mesh() const {
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    int    iV0   = mesh.face_vrtx(iF,0) ;
    int    iV1   = mesh.face_vrtx(iF,1) ;
    double xV[2] = { mesh.coords_V(iV0,0), mesh.coords_V(iV1,0) } ;
    double yV[2] = { mesh.coords_V(iV0,1), mesh.coords_V(iV1,1) } ;
    psd.draw_line( xV[0], yV[0], xV[1], yV[1] ) ;
  }
}
void mesh2D_psout::draw_median_dual_mesh() const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double x0 = mesh.coords_R(iR,0) ;
    double y0 = mesh.coords_R(iR,1) ;
    vector<int> flist ;
    mesh.get_regn_face( iR, flist ) ;
    int nRF = flist.size() ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = flist[ilF] ;
      double x1 = mesh.coords_F(iF,0) ;
      double y1 = mesh.coords_F(iF,1) ;
      psd.draw_line( x0, y0, x1, y1 ) ;
    }
  }
}
void mesh2D_psout::draw_center_dual_mesh() const {
  int nF = mesh.n_face () ;
  for ( int iF=0; iF<nF; ++iF ) {
    int    iR0 = mesh.face_regn(iF,0) ;
    double xR0 = mesh.coords_R(iR0,0) ;
    double yR0 = mesh.coords_R(iR0,1) ;    
    if ( mesh.is_boundary_face(iF) ) {
      double xbF = mesh.coords_F(iF,0) ; ;
      double ybF = mesh.coords_F(iF,1) ; ;
     psd.draw_line( xR0, yR0, xbF, ybF ) ;
    } else {
      int    iR1 = mesh.face_regn(iF,1) ;
      double xR1 = mesh.coords_R(iR1,0) ;
      double yR1 = mesh.coords_R(iR1,1) ;
      psd.draw_line( xR0, yR0, xR1, yR1 ) ;
    }
  }
}
void mesh2D_psout::draw_EE_graph () const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> flist ;
    mesh.get_regn_face(iR,flist) ;
    int nRF = flist.size() ;
    for ( int il0=0; il0<nRF; ++il0 ) {
      int    il1   = (il0+1)%nRF ;
      int    iF0   = flist[il0] ;
      int    iF1   = flist[il1] ;
      double xF[2] = { mesh.coords_F(iF0,0), mesh.coords_F(iF1,0) } ; 
      double yF[2] = { mesh.coords_F(iF0,1), mesh.coords_F(iF1,1) } ; 
      psd.draw_line( xF[0], yF[0], xF[1], yF[1] ) ;
    }
  }
}
void mesh2D_psout::draw_ET_graph() const {
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    double xFm = mesh.coords_F(iF,0) ;
    double yFm = mesh.coords_F(iF,1) ;
    int    iR0 = mesh.face_regn(iF,0) ;
    double xR0 = mesh.coords_R(iR0,0) ;
    double yR0 = mesh.coords_R(iR0,1) ;
    psd.draw_line( xFm, yFm, xR0, yR0 ) ;
    if ( mesh.is_internal_face(iF) ) {
      int    iR1 = mesh.face_regn(iF,1) ;
      double xR1 = mesh.coords_R(iR1,0) ;
      double yR1 = mesh.coords_R(iR1,1) ;
      psd.draw_line( xFm, yFm, xR1, yR1 ) ;
    }
  }
}
void mesh2D_psout::draw_TT_graph() const {
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double x0 = mesh.coords_R(iR,0) ;
    double y0 = mesh.coords_R(iR,1) ;
    vector<int> flist ;
    mesh.get_regn_face(iR,flist) ;
    int nRF = flist.size() ;
    for ( int il=0; il<nRF; ++il ) {
      int    iF = flist[il] ;
      double x1 = mesh.coords_F(iF,0) ;
      double y1 = mesh.coords_F(iF,1) ;
      psd.draw_line( x0, y0, x1, y1 ) ;
    }
  }
}
void mesh2D_psout::draw_bbox() const {
  double x0, y0, x1, y1 ;
  mesh.bbox( x0, y0, x1, y1 ) ;
  //psd.draw_rect( x0, y0, x1, y1 ) ;
  psd.draw_line( x0, y0, x1, y0 ) ;
  psd.draw_line( x1, y0, x1, y1 ) ;
  psd.draw_line( x1, y1, x0, y1 ) ;
  psd.draw_line( x0, y1, x0, y0 ) ;
}
void mesh2D_psout::draw_boundary() const {
  for ( int ilF=0; ilF<mesh.n_bface(); ++ilF ) {
    int ibF = mesh.get_bnd_face(ilF) ;
    int iV0 = mesh.face_vrtx(ibF,0) ;
    int iV1 = mesh.face_vrtx(ibF,1) ;
    double xV[2] = { mesh.coords_V(iV0,0), mesh.coords_V(iV1,0) } ;
    double yV[2] = { mesh.coords_V(iV0,1), mesh.coords_V(iV1,1) } ;
    psd.draw_line( xV[0], yV[0], xV[1], yV[1] ) ;
  }
}
void mesh2D_psout::draw_filled_regions( vector<bool> & regn_flag_list ) const {
  assert( regn_flag_list.size()==mesh.n_region() ) ;
  bool first_shot = bool(true) ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    if ( regn_flag_list[iR] ) {
      if ( first_shot ) {
	psd.dump_string( "% ----- filled regions" ) ;
	psd.set_rgb_colour( "LIGHTFILL" ) ;
	first_shot = bool(false) ;
      }
      double xR = mesh.coords_R(iR,0) ;
      double yR = mesh.coords_R(iR,1) ;
      vector<int> regn_vlist ;
      mesh.get_regn_vrtx( iR, regn_vlist ) ;
      int nRV = regn_vlist.size() ;
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	int iV0 = regn_vlist[ilV] ;
	int iV1 = regn_vlist[ (ilV+1)%nRV ] ;
	double xV0 = mesh.coords_V( iV0, 0 ) ;
	double yV0 = mesh.coords_V( iV0, 1 ) ;
	double xV1 = mesh.coords_V( iV1, 0 ) ;
	double yV1 = mesh.coords_V( iV1, 1 ) ;
	psd.fill_tria( xR, yR, xV0, yV0, xV1, yV1 );
      }
    }
  }
}
//draw_vrtx_pos
void mesh2D_psout::draw_vrtx_pos() const  {
  psd.dump_string( "% ----- vrtx numbers" ) ;
  psd.set_rgb_colour( "BLACK" ) ;
  psd.set_font( "Vrtx" ) ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    double xV = mesh.coords_V(iV,0) ;
    double yV = mesh.coords_V(iV,1) ;
    psd.draw_number( iV+get_offset(), xV, yV ) ;
  }
}
//draw_edge_pos
void mesh2D_psout::draw_edge_pos() const {
  psd.dump_string( "% ----- edge numbers" ) ;
  psd.set_rgb_colour( "RED" ) ;
  psd.set_font( "Edge" ) ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    double xF = mesh.coords_F(iF,0) ;
    double yF = mesh.coords_F(iF,1) ;
    psd.draw_number( iF+get_offset(), xF, yF ) ;
  }
}
//draw_regn_pos
void mesh2D_psout::draw_regn_pos() const {
  psd.dump_string( "% ----- regn numbers" ) ;
  psd.set_rgb_colour( "BLUE" ) ;
  psd.set_font( "Poly" ) ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    double xR = mesh.coords_R(iR,0) ;
    double yR = mesh.coords_R(iR,1) ;
    psd.draw_number( iR+get_offset(), xR, yR ) ;
  }
}
//draw_regn_bullet
void mesh2D_psout::draw_regn_bullet( int size, string col_str ) const {
  psd.dump_string( "% ----- regn bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iR=0; iR<mesh.n_region(); ++iR ) {
      double xR = mesh.coords_R(iR,0) ;
      double yR = mesh.coords_R(iR,1) ;
      // MODIFIED HMM position -------------------------
      if ( false ) {
	vector<int> R_vlist ;
	mesh.get_regn_vrtx( iR, R_vlist ) ;
	int nRV = R_vlist.size() ;
	double xRc = 0. ; 
	double yRc = 0. ;
	double sum = 0.; 
	for ( int ilV=0; ilV<3; ++ilV ) {
	  int iV = R_vlist[ilV] ;
	  xRc += mesh.coords_V( iV, 0 ) ;
	  yRc += mesh.coords_V( iV, 1 ) ;
	  sum += 1. ;
	}
	xRc /= sum ;
	yRc /= sum ;
	psd.draw_ellp( xRc, yRc, size ) ;
      } else {
	psd.draw_ellp( xR, yR, size ) ;
      }
    }
  }
}
//draw_vrtx_pos 
void mesh2D_psout::draw_int_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
      if ( mesh.is_internal_vrtx(iV) ) {
	double xV = mesh.coords_V(iV,0) ;
	double yV = mesh.coords_V(iV,1) ;
	psd.draw_ellp( xV, yV, size ) ;
      }
    }
  }
}   
void mesh2D_psout::draw_bnd_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
      if ( mesh.is_boundary_vrtx(iV) ) {
	double xV = mesh.coords_V(iV,0) ;
	double yV = mesh.coords_V(iV,1) ;
	psd.draw_ellp( xV, yV, size ) ;
      }
    }
  }
} 
void mesh2D_psout::draw_vrtx_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- vrtx bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
      double xV = mesh.coords_V(iV,0) ;
      double yV = mesh.coords_V(iV,1) ;
      psd.draw_ellp( xV, yV, size ) ;
    }
  }
}
//draw_vrtx_pos
void mesh2D_psout::draw_int_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iF=0; iF<mesh.n_face(); ++iF ) {
      if ( mesh.is_internal_face(iF) ) {
	double xF = mesh.coords_F(iF,0) ;
	double yF = mesh.coords_F(iF,1) ;
	psd.draw_ellp( xF, yF, size ) ;
      }
    }
  }
}
void mesh2D_psout::draw_bnd_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iF=0; iF<mesh.n_face(); ++iF ) {
      if ( mesh.is_boundary_face(iF) ) {
	double xF = mesh.coords_F(iF,0) ;
	double yF = mesh.coords_F(iF,1) ;
	psd.draw_ellp( xF, yF, size ) ;
      }
    }
  }
}
void mesh2D_psout::draw_face_bullet( int size, string col_str ) const  {
  psd.dump_string( "% ----- face bullets" ) ;
  psd.set_rgb_colour( col_str.c_str() ) ;
  psd.dump_string( "" ) ;
  if ( size>0 ) {
    for ( int iF=0; iF<mesh.n_face(); ++iF ) {
      double xF = mesh.coords_F(iF,0) ;
      double yF = mesh.coords_F(iF,1) ;
      psd.draw_ellp( xF, yF, size ) ;
    }
  }
}
void mesh2D_psout::draw_mesh_0( string mesh_label_str, bool draw_numbers ) const {
  setup_bb_mesh() ;                  // init bounding box

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ; 
  draw_primal_mesh() ;

  if ( draw_numbers ) {
    draw_vrtx_pos () ;              // draw vrtx/edge/regn numbers
    draw_edge_pos   () ;
    draw_regn_pos   () ;
  }

  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;
}
void mesh2D_psout::draw_mesh_1( string mesh_label_str, bool draw_numbers ) const {
  setup_bb_mesh() ;                 // init bounding box
  draw_bbox() ;                     // draw bounding box

  psd.set_rgb_colour( "BROWN" ) ;   // draw dual mesh
  psd.set_line_type ( "DASH1" ) ; 
  draw_median_dual_mesh() ;

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ;
  draw_primal_mesh() ;
  
  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
    draw_edge_pos() ;
    draw_regn_pos() ;
  }

  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;
}

void mesh2D_psout::draw_mesh( string mesh_label_str, bool draw_numbers ) const {
  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 750 TLS") ) ;

  setup_bb_mesh() ;                 // init bounding box
  //draw_bbox() ;                   // draw bounding box

#if 0
  psd.set_rgb_colour( "BLACK" ) ;   // draw dual mesh
  psd.set_line_type ( "DASH1" ) ;
  draw_center_dual_mesh() ;

  draw_regn_bullet    (3) ;
  draw_bnd_vrtx_bullet(3) ;
  draw_bnd_face_bullet(3) ;
#endif

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ;
  draw_primal_mesh() ;

  //draw_vrtx_bullet(2) ;

  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
    draw_edge_pos() ;
    draw_regn_pos() ;
    draw_vrtx_bullet(2) ;
    //draw_face_bullet(3) ;
  }
}
  
void mesh2D_psout::draw_mesh( string mesh_label_str, vector<bool> regn_flag_list,  
			      bool draw_numbers ) const {

  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 800 TLS") ) ;

  setup_bb_mesh() ;                 // init bounding box
  //draw_bbox() ;                   // draw bounding box

#if 0
  psd.set_rgb_colour( "BLACK" ) ;   // draw dual mesh
  psd.set_line_type ( "DASH1" ) ;
  draw_center_dual_mesh() ;

  draw_regn_bullet    (3) ;
  draw_bnd_vrtx_bullet(3) ;
  draw_bnd_face_bullet(3) ;
#endif

  draw_filled_regions( regn_flag_list ) ;

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ; 
  draw_primal_mesh() ;

  draw_vrtx_bullet(2) ;

  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
    draw_edge_pos() ;
    draw_regn_pos() ;
    draw_vrtx_bullet(2) ;
    draw_face_bullet(3) ;
  }
}

void mesh2D_psout::draw_meshless( string mesh_label_str, bool draw_numbers, bool draw_dualmesh ) const {
  psd.set_rgb_colour( "BLACK" ) ;   // set figure title
  psd.set_font( "Mesh" ) ;
  psd.dump_string( string("(") + mesh_label_str + string(") 50 750 TLS") ) ;

  setup_bb_mesh() ;                 // init bounding box
  draw_bbox() ;                   // draw bounding box

  if ( draw_dualmesh ) {
    psd.set_rgb_colour( "BLACK" ) ;   // draw dual mesh
    psd.set_line_type ( "DASH1" ) ;
    draw_center_dual_mesh() ;
    
    draw_regn_bullet    (3) ;
    draw_bnd_vrtx_bullet(3) ;
    draw_bnd_face_bullet(3) ;
  }

  psd.set_rgb_colour( "BLACK" ) ;   // draw primal mesh
  psd.set_line_type ( "SOLID" ) ;
  //draw_primal_mesh() ;

  draw_vrtx_bullet(3) ;

  if ( draw_numbers ) {
    draw_vrtx_pos() ;              // draw vrtx/edge/regn numbers
  }
}

# endif // end of _MESH_2D_PSOUT
