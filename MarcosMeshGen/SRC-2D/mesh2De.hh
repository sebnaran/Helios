#ifndef _MESH2DV_HH
#define _MESH2DV_HH

#include <limits>

namespace position_flags { 
  // position flags 
  const int UNSET = -9999999 ;
  const int flag_INTERNAL =    0 ;
  const int flag_BOUNDARY =   -1 ;
} ;

using namespace position_flags ;

class mesh_2Dv {

  // ---------------------------------------------
  // datasets
  // ---------------------------------------------
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
    int & operator() ( int i, int j ) {
      assert( 0<=i && i<nrow ) ;
      assert( 0<=j && j<iset[i].size() ) ;
      return iset[i][j] ;
    }
  } ;
  // ---------------------------------------------

  // friendships and type declarations  
  friend class mesh2Dv_builder ;
  friend class mesh2Dv_printer ;
  
protected:
  static const int DIM=2 ;
  static const int nFR=2 ;
  static const int nFV=2 ;

  int nV, nF, nR ;
  string mesh_name ;

  // coordinates
  vector<double> vrtx_coords ;
  vector<double> regn_coords ;

  // bounding box & mesh size
  bool   bb_status ;
  double bb_min[DIM], bb_max[DIM] ;
  
  bool   ms_status ;
  double ms_hmax, ms_hmin, ms_hsqr ;
  
  // datasets
  dataset RegnFace ;
  dataset VrtxFace ;

  // vectors
  vector<int> FaceVrtx ;
  vector<int> FaceRegn ;

  // lists of boundary items
  vector<int> bnd_vrtx ;
  vector<int> bnd_face ;
  vector<int> bnd_regn ;

  // external flags
  vector<int> fV ;
  vector<int> fF ;
  vector<int> fR ;

  // aux geometrical quantities
  vector<double> regn_area   ;
  vector<double> face_length ;

  // private method
  void shrink_list( vector<int> & tmp_list ) ;

  // reset status
  void reset_status_flags() {
    bb_status = false ;
    ms_status = false ;
  }

public:
  mesh_2Dv() : bb_status(false), ms_status(false), mesh_name("mesh-2D") {}
  ~mesh_2Dv() {}

  // ACCESS METHODS: cardinality of mesh item sets
  inline int n_region() { return nR ; }
  inline int n_face  () { return nF ; }
  inline int n_vertex() { return nV ; }

  inline int n_bregion() { return bnd_regn.size() ; }
  inline int n_bface  () { return bnd_face.size() ; }
  inline int n_bvertex() { return bnd_vrtx.size() ; }

  inline int get_bnd_regn( int ilR ) ;
  inline int get_bnd_face( int ilF ) ;
  inline int get_bnd_vrtx( int ilV ) ;

  // TOPOLOGICAL METHODS, (i_glob, i_loc)
  inline int regn_face( int iR, int ilF ) ;
  inline int face_vrtx( int iF, int ilV ) ;

  inline int face_regn( int iF, int ilR ) ;
  inline int vrtx_face( int iV, int ilF ) ;

  inline int vrtx_vrtx( int iV, int ilV ) ; // ISO vrtx_edge

  inline int n_regn_face( int iR ) ;
  inline int n_face_vrtx( int iF ) ;

  inline int n_face_regn( int iF ) ;
  inline int n_vrtx_face( int iV ) ;

  inline int n_vrtx_vrtx( int iV ) ;

  inline bool ok_regn_face( int iR, int ilF ) ;
  inline bool ok_vrtx_face( int iV, int ilF ) ;
  
  // LISTS OF MESH ITEMS, ( i_glob, list )
  // ...for regions
  inline void get_regn_regn( int iR, vector<int> & rlist ) ;
  inline void get_regn_face( int iR, vector<int> & flist ) ;
  inline void get_regn_vrtx( int iR, vector<int> & vlist ) ;

  // ...for faces
  inline void get_face_regn( int iF, vector<int> & rlist ) ;
  inline void get_face_face( int iF, vector<int> & flist ) ;
  inline void get_face_vrtx( int iF, vector<int> & vlist ) ;

  // ...for vertices
  inline void get_vrtx_regn( int iV, vector<int> & rlist ) ;
  inline void get_vrtx_face( int iV, vector<int> & flist ) ;
  inline void get_vrtx_vrtx( int iV, vector<int> & vlist ) ;
  
  // Logical Methods for detecting boundary items
  inline bool is_boundary_vrtx( int iV ) { return binary_search( bnd_vrtx.begin(), bnd_vrtx.end(), iV ) ; }
  inline bool is_boundary_face( int iF ) { return binary_search( bnd_face.begin(), bnd_face.end(), iF ) ; }
  inline bool is_boundary_regn( int iR ) { return binary_search( bnd_regn.begin(), bnd_regn.end(), iR ) ; }

  inline bool is_internal_vrtx( int iV ) { return !is_boundary_vrtx(iV) ; }
  inline bool is_internal_face( int iF ) { return !is_boundary_face(iF) ; }
  inline bool is_internal_regn( int iR ) { return !is_boundary_regn(iR) ; }

  // return local position of face iF inside regn iR or -1 if iF is not a face of iR
  int get_local_regn_face( int iR, int iF ) {
    assert( 0<=iR && iR<nR ) ;
    assert( 0<=iF && iF<nF ) ;
    int rval = -1 ;
    vector<int> regn_flist ;
    get_regn_face( iR, regn_flist ) ;
    int nRF = regn_flist.size() ;
    for ( int ilF=0; ilF<nRF && rval==-1; ++ilF ) {
      if ( regn_flist[ilF]==iF ) {
	rval = ilF ;
      }
    }
    return rval ;
  }

  // Geometrical Methods
  inline double get_nor( int iF, int s ) ;
  inline double get_tng( int iF, int s ) ;

  // Geometrical measures
  inline double get_regn_measure( int iR ) ;
  inline double get_face_measure( int iF ) ;

  // eval bounding box
  inline double eval_bbox( string retstr ) ;
  inline void   eval_bbox() ;

  inline double min_coords( int s ) { 
    assert( 0<=s && s<DIM ) ;
    if ( !bb_status ) { eval_bbox() ; }
    return bb_min[s] ;
  }
  inline double max_coords( int s ) { 
    assert( 0<=s && s<DIM ) ;
    if ( !bb_status ) { eval_bbox() ; }
    return bb_max[s] ;
  }

  inline void change_bbox( double new_xmin, double new_ymin, double new_xmax, double new_ymax ) ;
  inline void bbox( double & xmin, double & ymin, double & xmax, double & ymax ) ;

  // used by problems in pblm
  inline double xmin()  { return bb_status ? bb_min[0] : eval_bbox("xmin") ; }
  inline double ymin()  { return bb_status ? bb_min[1] : eval_bbox("ymin") ; }
  inline double xmax()  { return bb_status ? bb_max[0] : eval_bbox("xmax") ; }
  inline double ymax()  { return bb_status ? bb_max[1] : eval_bbox("ymax") ; }

  // Geometrical Methods (with some problems)
  inline double coords_V( int iV, int s ) ;
  inline double coords_F( int iF, int s ) ;
  inline double coords_R( int iR, int s ) ;

  // Mesh sizes: hmax/hmin/hsqr
  inline double h_max() ;
  inline double h_min() ;
  inline double h_sqr() ;
  inline double h_regn( int iR ) ;
  inline double h_regn_face( int iR ) ;
  inline double eval_mesh_size( string retstr="" ) ;

  // change coordinates of vertices
  inline void set_vrtx_coords( int iV, double xv, double yv ) {
    vrtx_coords[DIM*iV+0] = xv ;
    vrtx_coords[DIM*iV+1] = yv ;
  }

  // mesh name
  inline void   set_mesh_name( string _mesh_name ) { mesh_name = _mesh_name ; } 
  inline string get_mesh_name() { return mesh_name ; } 

  // reset external flags (useful for gmv option "explode")
  void set_fV( int iV, int new_fV ) { fV[iV] = new_fV ; }
  void set_fF( int iF, int new_fF ) { fF[iF] = new_fF ; }
  void set_fR( int iR, int new_fR ) { fR[iR] = new_fR ; }

  // get external flags
  int  get_fV( int iV ) { return fV[iV] ; }
  int  get_fF( int iF ) { return fF[iF] ; }
  int  get_fR( int iR ) { return fR[iR] ; }

  // --------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------
  //
  // EDGES HAVE NO IMPLEMENTATION, ONLY INTERFACE
  //
  // --------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------

  inline int n_edge  () { return nF ; }
  inline int n_bedge () { return bnd_face.size() ; }
  inline int get_bnd_edge( int ilE ) { return get_bnd_face(ilE) ; }

  // topological methods
  inline int regn_edge( int iR, int ilE ) { return regn_face(iR,ilE) ; }
  inline int edge_vrtx( int iE, int ilV ) { return face_vrtx(iE,ilV) ; }

  inline int edge_regn( int iE, int ilR ) { return face_regn(iE,ilR) ; }
  inline int vrtx_edge( int iV, int ilE ) { return vrtx_edge(iV,ilE) ; }

  inline int n_regn_edge( int iR ) { return n_regn_face(iR) ; }
  inline int n_edge_vrtx( int iE ) { return n_face_vrtx(iE) ; }

  inline int n_edge_regn( int iE ) { return n_face_regn(iE) ; }
  inline int n_vrtx_edge( int iV ) { return n_vrtx_edge(iV) ; }

  inline bool ok_regn_edge( int iR, int ilE ) { return ok_regn_face(iR,ilE) ; }
  inline bool ok_vrtx_edge( int iV, int ilE ) { return ok_vrtx_edge(iV,ilE) ; }
  
  // return list for regions
  inline void get_regn_edge( int iR, vector<int> & elist ) { get_regn_face(iR,elist) ; }

  // return lists for edges
  inline void get_edge_regn( int iE, vector<int> & rlist ) { get_edge_regn(iE,rlist) ; }
  inline void get_edge_edge( int iE, vector<int> & elist ) { get_edge_edge(iE,elist) ; }
  inline void get_edge_vrtx( int iE, vector<int> & vlist ) { get_edge_vrtx(iE,vlist) ; }

  // return list for vertices
  inline void get_vrtx_edge( int iV, vector<int> & elist ) { get_vrtx_face(iV,elist) ; }
  
  // Logical Methods for detecting boundary items
  inline bool is_boundary_edge( int iE ) { return is_boundary_face(iE) ; }
  inline bool is_internal_edge( int iE ) { return !is_boundary_face(iE) ; }

  int get_local_regn_edge( int iR, int iE ) {
    return get_local_regn_face( iR, iE ) ;
  }

  // misc
  inline double get_edge_measure( int iE ) { return get_face_measure(iE) ; }
  int  get_fE( int iE ) { return fF[iE] ; }

  // --------------------------------------------------------------------------------------
} ;
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
void mesh_2Dv :: shrink_list( vector<int> & tmp_list ) {
  sort( tmp_list.begin(), tmp_list.end() ) ;
  int k = 0 ;
  for ( int i=1; i<tmp_list.size(); ++i ) {
    if ( tmp_list[i]!=tmp_list[k] ) {
      tmp_list[++k]=tmp_list[i] ;
    }
  }
  tmp_list.resize(k+1) ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
  // TOPOLOGICAL METHODS:
int mesh_2Dv :: regn_face( int iR, int ilF ) {
  assert( 0<=iR  && iR<nR ) ;
  assert( 0<=ilF && ilF<RegnFace.size_loc(iR) ) ;
  return RegnFace(iR,ilF) ;
}
int mesh_2Dv :: face_vrtx( int iF, int ilV ) {
  assert( 0<=iF  && iF<nF  ) ;
  assert( ilV==0 || ilV==1 ) ;
  return FaceVrtx[nFV*iF+ilV] ;
}
int mesh_2Dv :: face_regn( int iF, int ilR ) {
  assert( 0<=iF  && iF<nF  ) ;
  assert( ilR==0 || ilR==1 ) ;
  return FaceRegn[nFR*iF+ilR] ;
}
int  mesh_2Dv :: vrtx_face( int iV, int ilF ) {
  assert( 0<=iV  && iV<nV ) ;
  assert( 0<=ilF && ilF<VrtxFace.size_loc(iV) ) ;
  return VrtxFace(iV,ilF) ;
}

// vrtx_vrtx ISO vrtx_face
int mesh_2Dv :: vrtx_vrtx( int iV, int ilF ) {
  assert( 0<=iV && iV<nV ) ;
  int iF  = VrtxFace(iV,ilF) ; 
  int iV0 = FaceVrtx[nFV*iF+0] ;
  int iV1 = FaceVrtx[nFV*iF+1] ;
  return iV0==iV ? iV1 : iV0 ;
}

int mesh_2Dv :: n_regn_face( int iR ) {
  assert( 0<=iR  && iR<nR ) ;
  return RegnFace.size_loc(iR) ;
}
int mesh_2Dv :: n_face_vrtx( int iF ) {
  assert( 0<=iF  && iF<nF ) ;
  return nFV ;
}
// vrtx_vrtx ISO vrtx_face
int mesh_2Dv :: n_vrtx_vrtx( int iV ) {
  assert( 0<=iV  && iV<nV ) ;
  return VrtxFace.size_loc(iV) ;
}

int mesh_2Dv :: n_face_regn( int iF ) {
  assert( 0<=iF  && iF<nF ) ;
  return nFR ;
}
int mesh_2Dv :: n_vrtx_face( int iV ) {
  assert( 0<=iV  && iV<nV ) ;
  return VrtxFace.size_loc(iV) ;
}

// ok_regn_face == TRUE if face is pointing out of region
bool mesh_2Dv :: ok_regn_face( int iR, int ilF ) {
  assert( 0<=iR  && iR<nR ) ;
  assert( 0<=ilF && ilF<RegnFace.size_loc(iR) ) ;
  int iF = regn_face(iR,ilF) ;
  return iR==face_regn(iF,0) ;
}
// ok_vrtx_face == TRUE if face is pointing out of vertex
bool mesh_2Dv :: ok_vrtx_face( int iV, int ilF ) {
  assert( 0<=iV  && iV<nV ) ;
  assert( 0<=ilF && ilF<VrtxFace.size_loc(iV) ) ;
  int iF = vrtx_face(iV,ilF) ;
  return iV==face_vrtx(iF,0) ;  
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of regions
// the "UNSET" choice for boundary face preserves the order consistency
// with face and vertex lists
void mesh_2Dv :: get_regn_regn( int iR, vector<int> & rlist ) {
  assert( 0<=iR && iR<nR ) ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { 
    int iF  = RegnFace(iR,ilF) ;
    int iR0 = FaceRegn[nFR*iF+0] ;
    int iR1 = FaceRegn[nFR*iF+1] ;
    if      ( iR==iR1               ) { rlist.push_back( iR0 ) ; }
    else if ( iR==iR0 && iR1!=UNSET ) { rlist.push_back( iR1 ) ; }
    else if ( iR==iR0 && iR1==UNSET ) { rlist.push_back( iR1 ) ; } // boundary face --> UNSET
    else { // it should never happen
      MSG("-->>get_regn_regn: consistency error! "<<endl<<flush) ; 
      assert(false) ; 
    }
  }
}
void mesh_2Dv :: get_regn_face( int iR, vector<int> & flist ) {
  assert( 0<=iR && iR<nR ) ;
  flist.resize( RegnFace.size_loc(iR) ) ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { flist[ilF] = RegnFace(iR,ilF) ; }
}
void mesh_2Dv :: get_regn_vrtx( int iR, vector<int> & vlist ) {
  assert( 0<=iR && iR<nR ) ;
  vlist.resize( 0 ) ; 
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { 
    int iF = RegnFace(iR,ilF) ;
    if ( ok_regn_face(iR,ilF) ) { vlist.push_back( FaceVrtx[nFV*iF+0] ) ; }
    else                        { vlist.push_back( FaceVrtx[nFV*iF+1] ) ; }
  }
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of faces
void mesh_2Dv :: get_face_regn( int iF, vector<int> & rlist ) {
  assert( 0<=iF && iF<nF ) ;
  rlist.resize(0) ;
  rlist.push_back(FaceRegn[nFR*iF+0]) ;
  if ( FaceRegn[nFR*iF+1]!=UNSET ) {
    rlist.push_back(FaceRegn[nFR*iF+1]) ;
  }
}
void mesh_2Dv :: get_face_face( int iF, vector<int> & flist ) {
  assert( 0<=iF  && iF<nF ) ;
  flist.resize( 0 ) ;
  for ( int ilV=0; ilV<nFV; ++ilV ) {
    int iV = FaceVrtx[nFV*iF+ilV] ;
    for ( int ilF=0; ilF<VrtxFace.size_loc(iV); ++ilF ) {
      int jF = VrtxFace(iV,ilF) ;
      if ( iF != jF ) { flist.push_back( jF ) ; }
    }
  } 
}
void mesh_2Dv :: get_face_vrtx( int iF, vector<int> & vlist ) {
  assert( 0<=iF && iF<nF ) ;
  vlist.resize( 2 ) ;
  vlist[0] = FaceVrtx[nFV*iF+0] ;
  vlist[1] = FaceVrtx[nFV*iF+1] ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of vertices
void mesh_2Dv :: get_vrtx_regn( int iV, vector<int> & rlist ) {
  assert( 0<=iV && iV<nV ) ;
  vector<int> tmp_vec ;  
  for ( int ilF=0; ilF<VrtxFace.size_loc(iV); ++ilF ) { 
    int iF = VrtxFace(iV,ilF) ; 
    tmp_vec.push_back( FaceRegn[nFR*iF+0] ) ;
    if ( FaceRegn[nFR*iF+1]!=UNSET ) {
      tmp_vec.push_back( FaceRegn[nFR*iF+1] ) ;
    }
  }
  shrink_list( tmp_vec ) ;
  rlist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<rlist.size(); ++i ) { rlist[i]=tmp_vec[i] ; }
}
void mesh_2Dv :: get_vrtx_face( int iV, vector<int> & flist ) {
  assert( 0<=iV && iV<nV ) ;
  flist.resize( 0 ) ;
  for ( int ilF=0; ilF<VrtxFace.size_loc(iV); ++ilF ) { 
    flist.push_back( VrtxFace(iV,ilF) ) ; 
  }
}
void mesh_2Dv :: get_vrtx_vrtx( int iV, vector<int> & vlist ) {
  assert( 0<=iV && iV<nV ) ;
  vlist.resize( VrtxFace.size_loc(iV) ) ;
  for ( int ilF=0; ilF<VrtxFace.size_loc(iV); ++ilF ) { 
    int iF  = VrtxFace(iV,ilF) ; 
    int iV0 = FaceVrtx[nFV*iF+0] ;
    int iV1 = FaceVrtx[nFV*iF+1] ;
    vlist[ilF] = iV0==iV ? iV1 : iV0 ;
  }
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// (new)
double mesh_2Dv :: get_regn_measure( int iR ) {
  assert( 0<=iR && iR<nR ) ;
  assert( regn_area.size()==nR ) ;
  return regn_area[iR] ;
}
double mesh_2Dv :: get_face_measure( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  assert( face_length.size()==nF ) ;
  return face_length[iF] ;
}
double mesh_2Dv :: get_nor( int iF, int s ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=s  && s<DIM ) ;
  int    iV0   = face_vrtx(iF,0) ;
  int    iV1   = face_vrtx(iF,1) ;
  double len_F = get_face_measure(iF) ;
  double sgn_nor[DIM] = { +1., -1. } ;
  return sgn_nor[s] * ( coords_V( iV1, (s+1)%DIM )-coords_V( iV0, (s+1)%DIM ) )/len_F ;
}
double mesh_2Dv :: get_tng( int iF, int s ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=s  && s<DIM ) ;
  int    iV0   = face_vrtx(iF,0) ;
  int    iV1   = face_vrtx(iF,1) ;
  double len_F = get_face_measure(iF) ;
  return ( coords_V(iV1,s)-coords_V(iV0,s) )/len_F ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
int mesh_2Dv :: get_bnd_regn( int ilR ) { 
  assert( 0<=ilR && ilR<bnd_regn.size() ) ;
  return bnd_regn[ilR] ;
}
int mesh_2Dv :: get_bnd_face( int ilF ) { 
  assert( 0<=ilF && ilF<bnd_face.size() ) ;
  return bnd_face[ilF] ;
}
int mesh_2Dv :: get_bnd_vrtx( int ilV ) { 
  assert( 0<=ilV && ilV<bnd_vrtx.size() ) ;
  return bnd_vrtx[ilV] ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
double mesh_2Dv :: coords_R( int iR, int k ) {
  assert( 0<=iR && iR<nR ) ;
  assert( 0<=k  && k<DIM ) ;
  return regn_coords[ DIM*iR+k ] ;
}
double mesh_2Dv :: coords_V( int iV, int k ) {
  assert( 0<=iV && iV<nV ) ;
  assert( 0<=k  && k<DIM ) ;
  return vrtx_coords[DIM*iV+k] ;
}
double mesh_2Dv :: coords_F( int iF, int k ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=k  && k<DIM ) ;
  double retval = 0. ; 
  for ( int ilV=0; ilV<nFV; ++ilV ) { 
    int iV = FaceVrtx[nFV*iF+ilV] ;
    retval += coords_V(iV,k) ;
  }
  return retval/double(nFV) ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void mesh_2Dv :: eval_bbox() {
  bb_status = true ;
  for ( int s=0; s<DIM; ++s ) { 
    bb_min[s] = +1e+99 ;
    bb_max[s] = -1e+99 ;
  }
  for ( int iV=0; iV<n_vertex(); ++iV ) {
    for ( int s=0; s<DIM; ++s ) {
      bb_min[s] = min( bb_min[s], vrtx_coords[ DIM*iV+s ] ) ;
      bb_max[s] = max( bb_max[s], vrtx_coords[ DIM*iV+s ] ) ;
    }
  }
}
double mesh_2Dv :: eval_bbox( string retstr ) {
  bb_status = true ;
  for ( int s=0; s<DIM; ++s ) { 
    bb_min[s] = +1e+99 ;
    bb_max[s] = -1e+99 ;
  }
  for ( int iV=0; iV<n_vertex(); ++iV ) {
    for ( int s=0; s<DIM; ++s ) {
      bb_min[s] = min( bb_min[s], vrtx_coords[ DIM*iV+s ] ) ;
      bb_max[s] = max( bb_max[s], vrtx_coords[ DIM*iV+s ] ) ;
    }
  }
  double retval ;
  if      ( retstr=="xmin" ) { retval = bb_min[0] ; }
  else if ( retstr=="ymin" ) { retval = bb_min[1] ; }
  else if ( retstr=="xmax" ) { retval = bb_max[0] ; }
  else if ( retstr=="ymax" ) { retval = bb_max[1] ; }
  else                       { retval = 0.        ; }
  return retval ;
}

void mesh_2Dv :: bbox( double & xmin, double & ymin, double & xmax, double & ymax ) {
  if ( !bb_status ) { eval_bbox() ; }
  xmin = bb_min[0] ;
  ymin = bb_min[1] ;
  xmax = bb_max[0] ;
  ymax = bb_max[1] ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
double mesh_2Dv :: eval_mesh_size( string retstr ) {
  ms_status = true ;
  ms_hsqr = 0. ;
  for ( int iR=0; iR<n_region(); ++iR ) {
    ms_hsqr = max( ms_hsqr, sqrt( regn_area[iR] ) ) ;
  }
  ms_hmin = 1.e+99 ;
  for ( int iF=0; iF<n_face(); ++iF ) {
    ms_hmin = min( ms_hmin, face_length[iF] ) ;
  }

  // evaluate max distance of a polygon (J. Droniou's suggestion)
  ms_hmax = 0. ;
  for ( int iR=0; iR<n_region(); ++iR ) {
    double max_dist_R = 0. ;
    vector<int> R_vlist ;
    get_regn_vrtx( iR, R_vlist ) ;
    int nRV = R_vlist.size() ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int    iV  = R_vlist[ilV] ;
      double xiV = coords_V( iV, 0 ) ;
      double yiV = coords_V( iV, 1 ) ;
      for ( int jlV=ilV+1; jlV<nRV; ++jlV ) {
	int    jV  = R_vlist[jlV] ;
	double xjV = coords_V( jV, 0 ) ;
	double yjV = coords_V( jV, 1 ) ;
	double dist = sqrt( pow( xiV-xjV, 2 ) + pow( yiV-yjV, 2 ) ) ;
	max_dist_R = max( max_dist_R,  dist ) ;
      }
    }
    ms_hmax = max( ms_hmax,  max_dist_R ) ;
  }
  
  double retval ;
  if      ( retstr=="hmin" ) { retval = ms_hmin ; }
  else if ( retstr=="hmax" ) { retval = ms_hmax ; }
  else if ( retstr=="hsqr" ) { retval = ms_hsqr ; }
  else                       { retval = 0.      ; }
  return retval ;
}
inline double mesh_2Dv :: h_sqr() { return ms_status ? ms_hsqr : eval_mesh_size("hsqr") ; }
inline double mesh_2Dv :: h_max() { return ms_status ? ms_hmax : eval_mesh_size("hmax") ; }
inline double mesh_2Dv :: h_min() { return ms_status ? ms_hmin : eval_mesh_size("hmin") ; }

double mesh_2Dv :: h_regn_face( int iR ) {
  double retval = 0. ;
  vector<int> R_flist ;
  get_regn_face( iR, R_flist ) ;
  int nRF = R_flist.size() ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = R_flist[ilF] ;
    retval = max( retval, face_length[iF] ) ;
  }
  return retval ;
}

double mesh_2Dv :: h_regn( int iR ) {
  // evaluare max distance of a polygon (from a Jerome suggestion)
  double retval = 0. ;
  vector<int> R_vlist ;
  get_regn_vrtx( iR, R_vlist ) ;
  int nRV = R_vlist.size() ;
  for ( int ilV=0; ilV<nRV; ++ilV ) {
    int    iV  = R_vlist[ilV] ;
    double xiV = coords_V( iV, 0 ) ;
    double yiV = coords_V( iV, 1 ) ;
    for ( int jlV=ilV+1; jlV<nRV; ++jlV ) {
      int    jV  = R_vlist[jlV] ;
      double xjV = coords_V( jV, 0 ) ;
      double yjV = coords_V( jV, 1 ) ;
      double dist = sqrt( pow( xiV-xjV, 2 ) + pow( yiV-yjV, 2 ) ) ;
      retval = max( retval,  dist ) ;
    }
  }
  return retval ;
}

void mesh_2Dv :: change_bbox( double new_xmin, double new_ymin, double new_xmax, double new_ymax ) {
  assert( new_xmax>new_xmin ) ;
  assert( new_ymax>new_ymin ) ;
  
  if ( !bb_status ) { eval_bbox() ; }

  double bb_xmin = bb_min[0] ;
  double bb_ymin = bb_min[1] ;
  double bb_xmax = bb_max[0] ;
  double bb_ymax = bb_max[1] ;

  double xden = bb_xmax-bb_xmin ;
  double yden = bb_ymax-bb_ymin ;
  double ascl = (new_xmax-new_xmin)*(new_ymax-new_ymin)/( xden*yden ) ;

  for ( int iV=0; iV<nV; ++iV ) {
    double xV = vrtx_coords[ DIM*iV+0 ] ;
    double yV = vrtx_coords[ DIM*iV+1 ] ;
    vrtx_coords[ DIM*iV+0 ] = ( ( xV-bb_xmin ) * new_xmax + ( bb_xmax-xV ) * new_xmin ) / xden ;
    vrtx_coords[ DIM*iV+1 ] = ( ( yV-bb_ymin ) * new_ymax + ( bb_ymax-yV ) * new_ymin ) / yden ;
  }

  for ( int iR=0; iR<nR; ++iR ) {
    double xR = regn_coords[ DIM*iR+0 ] ;
    double yR = regn_coords[ DIM*iR+1 ] ;
    regn_coords[ DIM*iR+0 ] = ( ( xR-bb_xmin ) * new_xmax + ( bb_xmax-xR ) * new_xmin ) / xden ;
    regn_coords[ DIM*iR+1 ] = ( ( yR-bb_ymin ) * new_ymax + ( bb_ymax-yR ) * new_ymin ) / yden ;
    regn_area[iR] *= ascl ;
  }

  for ( int iF=0; iF<nF; ++iF ) {
    int    iV0 = face_vrtx(iF,0) ;
    int    iV1 = face_vrtx(iF,1) ;
    double xV0 = vrtx_coords[ DIM*iV0+0 ] ;
    double yV0 = vrtx_coords[ DIM*iV0+1 ] ;
    double xV1 = vrtx_coords[ DIM*iV1+0 ] ;
    double yV1 = vrtx_coords[ DIM*iV1+1 ] ;
    face_length[iF] = sqrt( pow(xV0-xV1,2)+pow(yV0-yV1,2) ) ;
  }
  bb_min[0] = new_xmin ;
  bb_min[1] = new_ymin ;
  bb_max[0] = new_xmax ;
  bb_max[1] = new_ymax ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

#endif // end of _MESH2DV_HH
