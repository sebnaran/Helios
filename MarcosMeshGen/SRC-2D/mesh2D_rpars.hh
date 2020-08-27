# ifndef _MESH_2D_RPARS_HH
# define _MESH_2D_RPARS_HH

//=================== standard includes

# include<iostream>
# include<fstream>
# include<string>
# include<vector>

using namespace std ;

// (using macros from my_macro.hh)

//=================== class Read_Mesh_Pars

class Read_Mesh2D_Pars {
protected:
  // flags
  int mesh_flag ;
  
  // mesh pars
  string mesh_name ;
  int    mesh_nx, mesh_ny, mesh_nlev ;
  double mesh_wx, mesh_wy ;
  
  // output pars
  int psgd_flag, numb_flag ;
  int prtr_flag ; 
  int outp_flag ; 

  // refinement par
  int    refn_flag  ;
  double refn_theta ;
  int    dual_flag, voro_flag ;

  // maximum  number of refinements
  int refn_max ;
  
  // ------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------
  // PRIVATE/PROTECTED METHODS
  // ------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------
  void set_defaults() {
    // flags
    mesh_flag = 0 ;
    
    // mesh pars
    mesh_name = "(unset)" ;
    mesh_nx = mesh_ny = mesh_nlev = 1 ;
    mesh_wx = mesh_wy = 0. ;
    
    // output parameters
    psgd_flag = -1 ; 
    numb_flag =  0 ; 
    prtr_flag = -1 ; 
    outp_flag = -1 ; 
    
    // splitting && refinement pars
    voro_flag  = 0 ;
    dual_flag  = 0 ;
    refn_flag  = 0 ;
    refn_theta = 0. ;
    refn_max   = 0 ;
  }
  bool read_line( string keywd, ifstream & inpf ) {
    bool retval = true ;
    if      ( keywd == "mesh" ) { read_mesh(inpf) ; }
    else if ( keywd == "psgd" ) { read_psgd(inpf) ; }
    else if ( keywd == "prtr" ) { read_prtr(inpf) ; }
    else if ( keywd == "outp" ) { read_outp(inpf) ; }
    else if ( keywd == "refn" ) { read_refn(inpf) ; }
    else if ( keywd == "mxrf" ) { read_mxrf(inpf) ; }
    else if ( keywd == "dual" ) { read_dual(inpf) ; }
    else                        { retval = false  ; }
    return retval ;
  }
  void read_mesh( ifstream & inpf ) {
    inpf >> mesh_flag ;
    int flag = abs(mesh_flag) ;
    if ( flag<100 ) { // input from file
      inpf >> mesh_name ;
    } else {
      inpf >> mesh_nx >> mesh_ny >> mesh_wx >> mesh_wy ;
    }
    inpf >> mesh_nlev ;
  }
  void read_psgd( ifstream & inpf ) { inpf >> psgd_flag >> numb_flag  ; }
  void read_prtr( ifstream & inpf ) { inpf >> prtr_flag               ; }
  void read_outp( ifstream & inpf ) { inpf >> outp_flag               ; }
  void read_refn( ifstream & inpf ) { inpf >> refn_flag >> refn_theta ; }
  void read_mxrf( ifstream & inpf ) { inpf >> refn_max                ; }
  void read_dual( ifstream & inpf ) { inpf >> dual_flag >> voro_flag  ; }
  
public:
  Read_Mesh2D_Pars () { set_defaults() ; }
  ~Read_Mesh2D_Pars() {}
  
  // access method for mesh pars
  int    get_mesh_flag () { return mesh_flag ; }
  int    get_mesh_nlev () { return mesh_nlev ; }
  int    get_mesh_nx   () { return mesh_nx   ; }
  int    get_mesh_ny   () { return mesh_ny   ; }
  double get_mesh_wx   () { return mesh_wx   ; }
  double get_mesh_wy   () { return mesh_wy   ; }
  string get_mesh_name () { return mesh_name ; }
  void   set_mesh_name ( string & _mesh_name ) { mesh_name = _mesh_name ; }
  
  // access methods for output pars
  int    get_psgd_flag () { return psgd_flag ; }
  int    get_numb_flag () { return numb_flag ; }
  int    get_prtr_flag () { return prtr_flag ; }
  int    get_outp_flag () { return outp_flag ; }
  
  // access methods for splitting && refinement pars
  int    get_max_refn  () { return refn_max  ; }
  int    get_refn_flag () { return refn_flag ; }
  double get_refn_theta() { return refn_theta; }
  int    get_dual_flag () { return dual_flag ; }
  int    get_voro_flag () { return voro_flag ; }
  
  void print_pars() {
    LINE(---) ;
    MSG("Read_Mesh_Pars:" << endl) ;
    LINE(---) ;
    PRT( get_mesh_flag() ) ;
    PRT( get_mesh_name() ) ;
    PRT( get_mesh_nx  () ) ;
    PRT( get_mesh_ny  () ) ;
    PRT( get_mesh_wx  () ) ;
    PRT( get_mesh_wy  () ) ;
    PRT( get_mesh_nlev() ) ;
    
    LINE(---) ;
    PRT( get_psgd_flag() ) ;
    PRT( get_numb_flag() ) ;
    PRT( get_prtr_flag() ) ;
    PRT( get_outp_flag() ) ;
    
    LINE(---) ;
    PRT( get_max_refn  () ) ;
    PRT( get_refn_flag () ) ;
    PRT( get_refn_theta() ) ;
    PRT( get_dual_flag () ) ; 
    PRT( get_voro_flag () ) ;
    LINE(---) ;
  }
} ;

# endif // end of _MESH_2D_RPARS_HH
