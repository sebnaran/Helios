# ifndef _READ_PARAMETERS_HH
# define _READ_PARAMETERS_HH

//=================== standard includes

# include<iostream>
# include<fstream>
# include<string>

using namespace std ;

//=================== includes base classes forming Read_Parameters

# include "../SRC-MISC/mesh2D_rpars.hh"

//=================== class Read_Parameters

class Read_Parameters : public Read_Mesh2D_Pars {

private:
  static istream & get_comments( istream & is ) {
    bool newline = true ;
    while( newline ) { 
      newline = false ;
      while( !is.eof() && (is.peek()==' ' || is.peek()=='\t' || is.peek()=='\n') ) { is.get() ; }
      while( !is.eof() && is.peek()=='#' ) {
	while( !is.eof() && is.peek()!='\n' ) { is.get() ; }
	newline = true ;
      }
    }
    return is ;
  }
  string fname ;
  void fatal_error( string fname ) {
    cerr << "fatal error in opening file " << fname << endl << flush ;
    assert(false) ;
  }
  void read_line( ifstream & inpf ) {
    string keywd = "" ;
    inpf >> get_comments >> keywd ;
    Read_Mesh2D_Pars::read_line( keywd, inpf ) ;
  }
public:
  Read_Parameters( string _fname ) : fname(_fname) {}
  ~Read_Parameters() {}
  void print_parameters() {
    Read_Mesh2D_Pars::print_pars() ;
  }
  void read_file( string _fname="" ) {
    string inp_fname = _fname=="" ? fname : _fname ;
    ifstream inpf( inp_fname.c_str() ) ;
    if ( !inpf.good() ) { fatal_error(inp_fname) ; }
    while ( !inpf.eof() ) { read_line( inpf ) ; }
    inpf.close()  ;
  }
} ;

// very useful shortcut (used everywhere in my application codes)
typedef Read_Parameters RPars ; 

# endif // end of _READ_PARAMETERS_HH

