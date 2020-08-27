#ifndef _DATASET_HH
#define _DATASET_HH

//#include <iostream>
//#include <cmath>
//#include <cassert>

//-----------------------------------------------------------------------------
// similar to full matrix data type
//-----------------------------------------------------------------------------

class dataset {
private:
  int nrow ;
  vector<int> * iset ;

public:
  dataset() {}
  ~dataset() {}
  void setup( int _nrow ) {
    nrow = _nrow ; 
    iset = new vector<int>[nrow] ;
  }
  void setup( int i, int nnz ) {
    assert( 0<=i && i<nrow ) ;
    iset[i].resize(nnz) ;
  }
  int size() const { 
    return nrow ; 
  }
  int size_loc( int i ) const { 
    assert( 0<=i && i<nrow ) ;
    return iset[i].size() ; 
  }
  //only for accessing indices
  int operator() ( int i, int j ) const {
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    return iset[i][j] ;
  }
  // "set"-methods
  void set_index( int i, int j, int ival ) {
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    iset[i][j] = ival ;
  }
} ;

#endif // end of _DATASET_HH
