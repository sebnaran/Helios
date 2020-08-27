#ifndef _MESH_2DV_BUILDER_HH
#define _MESH_2DV_BUILDER_HH

#include <iostream>

class mesh2Dv_builder {

  static const int DIM = 2 ;
  static const int nFV = 2 ;
  
  // DEBUG
  void print_vec_regn_face() ;

private:
  struct aux_struct {
  public:
    int  iGlb ;
    int  iloc ;
    aux_struct( int _iGlb, int _iloc ) : 
      iGlb(_iGlb), iloc(_iloc) {}
    ~aux_struct() {}
    bool operator< ( const aux_struct & S1 ) const { return iGlb<S1.iGlb ; }
  } ;
  
  typedef struct aux_struct regn_face_struct ;

private:
  mesh_2Dv & mesh ;

private:
  bool verbose_flag ;

private:
  // build aux data structure
  void build_face_vlist( vector< vector<regn_face_struct> > & vec_regn_face,
			 vector<int>                        & face_vlist,
			 vector<int>                        & regn_flist ) ;

  void build_regn_flist( vector< vector<regn_face_struct> > & vec_regn_face,     // output (construction)
			 vector<int>                        & face_vrtx_list,    // input
			 vector<int>                        & face_regn_list ) ; // input

  // mesh coordinates
  void set_coords( vector<double> & xV, vector<double> & yV ) ;
  void set_flags ( vector<int>    & fV, vector<int>    & fR ) ;

  // build primary data structures
  void build_RegnFace( vector< vector<regn_face_struct> > & vec_regn_face, int _nR ) ;
  void build_FaceVrtx( vector<int> & face_vlist ) ;

  // build transposed data structures
  void build_VrtxFace() ;
  void build_FaceRegn() ;

  // build boundary lists
  void build_boundary_lists() ;
  void shrink_list( vector<int> & tmp_list ) ;

  // build geometric quantities
  void setup_geom_factors   () ;
  void set_regn_geom_factors() ;
  void set_face_geom_factors() ;

  // additional 
  void build_RegnFace( vector<int> & regn_face_list, vector<int> & face_regn_list ) ;

public:
  mesh2Dv_builder( mesh_2Dv & _mesh ) : mesh(_mesh), verbose_flag(true) {}
  ~mesh2Dv_builder() {}

  void change_bbox( double new_xmin, double new_ymin, double new_xmax, double new_ymax ) ;
  void recompute_geom_factors( bool _verbose_flag=true ) { 
    verbose_flag = _verbose_flag ;
    setup_geom_factors() ; 
  }

  void reorder_faces( vector<int> & face_perm ) ;

  void build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, 
		       vector<int> & regn_vlist, vector<int> & fR ) ;

  void build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, 
		       vector<int> & face_vrtx_list, vector<int> & face_regn_list, vector<int> & fF ) ;

  void build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV,
		       vector<int> & regn_vrtx_list, vector<int> & regn_face_list,
		       vector<int> & face_vrtx_list, vector<int> & face_regn_list, 
		       vector<int> & fF ) ;
} ; 

void mesh2Dv_builder ::
build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, 
		vector<int> & regn_vlist, vector<int> & fR ) {
  DBGF(begin  ---> mesh2Dv_builder::build_the_mesh using regn_vlist <----) ;
  
  // -- set coordinates in mesh
  set_coords( xV, yV ) ;
  
  // -- core for building primary dataset
  int nR = fR.size() ;
  vector< vector<regn_face_struct> > vec_regn_face(nR) ;
  vector<int> face_vlist ;
  build_face_vlist( vec_regn_face, face_vlist, regn_vlist ) ;
  build_RegnFace( vec_regn_face, nR ) ;
  
  // set face-vertex data structure
  build_FaceVrtx( face_vlist ) ;
  
  // -- set external flags for regions and vertices
  set_flags( fV, fR ) ;
  
  // set face flags to UNSET
  int nF = face_vlist.back() ;
  mesh.fF.resize(nF) ;
  for ( int iF=0; iF<nF; ++iF ) { mesh.fF[iF] = UNSET ; }
  
  // -- tranpose datasets
  build_VrtxFace() ;
  build_FaceRegn() ;
  
  // -- build boundary lists
  build_boundary_lists() ;
  
  // -- compute/set last geometric quantities
  setup_geom_factors() ;
  
  // final reset of logical flags
  mesh.reset_status_flags() ;
  DBGF(end of ---> mesh2Dv_builder::build_the_mesh using regn_vlist <----) ;
}

// setup coordinates
void mesh2Dv_builder ::
set_coords( vector<double> & xV, vector<double> & yV ) {
  assert( xV.size()==yV.size() ) ;
  int nV = xV.size() ;
  mesh.nV = nV ;
  mesh.vrtx_coords.resize(nV*DIM) ;
  for ( int iV=0; iV<nV; ++iV ) {
    mesh.vrtx_coords[DIM*iV+0] = xV[iV] ;
    mesh.vrtx_coords[DIM*iV+1] = yV[iV] ;
  }
}
// setup external flags
void mesh2Dv_builder ::
set_flags ( vector<int> & fV, vector<int> & fR ) {
  // --
  int nV = fV.size() ;
  mesh.fV.resize(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { mesh.fV[iV] = fV[iV] ; }
  // --
  int nR = fR.size() ;
  mesh.fR.resize(nR) ;
  for ( int iR=0; iR<nR; ++iR ) { mesh.fR[iR] = fR[iR] ; }
  // --
}
// build up primary datasets
void mesh2Dv_builder :: build_RegnFace( vector< vector<regn_face_struct> > & vec_regn_face, int _nR ) {
  std::cout << "begin build_RegnFace" << std::endl << std::flush ;
  int nR = _nR ;
  mesh.nR = nR ;
  mesh.RegnFace.setup( nR ) ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = vec_regn_face[iR].size() ;
    mesh.RegnFace.setup(iR,nRF) ;
    for ( int il=0; il<nRF; ++il ) {
      int  iF  = vec_regn_face[iR][il].iGlb ;
      int  ilF = vec_regn_face[iR][il].iloc ;
      mesh.RegnFace(iR,ilF) = iF ;
    }
  }
  std::cout << "end-->build_RegnFace" << std::endl << std::flush ;
}
void mesh2Dv_builder :: build_FaceVrtx( vector<int> & face_vlist ) {
  std::cout << "begin build_FaceVrtx" << std::endl << std::flush ;
  const int nFV = 2 ;
  const int nF  = face_vlist.back() ;
  mesh.nF = nF ;
  mesh.FaceVrtx.resize( 2*nF ) ;
  int k = 0 ;
  for ( int iF=0; iF<nF; ++iF ) {
    int kV0 = face_vlist[k++] ;
    int kV1 = face_vlist[k++] ;
    mesh.FaceVrtx[nFV*iF+0]=kV0 ;
    mesh.FaceVrtx[nFV*iF+1]=kV1 ;
  }
  std::cout << "end-->build_FaceVrtx" << std::endl << std::flush ;
}
// build up transposed datasets
// two-steps implementations are needed because
// one needs the local size of each list of transposed items
void mesh2Dv_builder :: build_VrtxFace() {
  std::cout << "begin build_VrtxFace" << std::endl << std::flush ;
  int nV  = mesh.nV ;
  int nF  = mesh.nF ;
  int nFV = 2 ;
  vector< vector<aux_struct> > vec_aux(nV) ;
  for ( int iF=0; iF<nF; ++iF ) {
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      int iV = mesh.FaceVrtx[nFV*iF+ilV] ;
      assert( 0<=iV && iV<nV ) ;
      vec_aux[iV].push_back( aux_struct(iF,ilV) ) ;
    }
  }
  mesh.VrtxFace.setup( nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    int nVF = vec_aux[iV].size() ;
    mesh.VrtxFace.setup( iV, nVF ) ;
    for ( int ilF=0; ilF<nVF; ++ilF ) {
      int  iF   = vec_aux[iV][ilF].iGlb ;
      assert( 0<=iF && iF<nF ) ;
      mesh.VrtxFace(iV,ilF) = iF ;
    }
  }
  std::cout << "end-->build_VrtxFace" << std::endl << std::flush ;
}
void mesh2Dv_builder :: build_FaceRegn() {
  std::cout << "begin build_FaceRegn" << std::endl << std::flush ;
  int nF = mesh.nF ;
  int nR = mesh.nR ;
  vector< vector<aux_struct> > vec_aux(nF) ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = mesh.RegnFace.size_loc(iR) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = mesh.RegnFace(iR,ilF) ;
      vec_aux[iF].push_back( aux_struct(iR,ilF) ) ;
    }
  }
  mesh.FaceRegn.resize( nF*2, UNSET ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int nFR = vec_aux[iF].size() ;
    assert( nFR==1 || nFR==2 ) ;
    for ( int ilR=0; ilR<nFR; ++ilR ) {
      int iR = vec_aux[iF][ilR].iGlb ;
      mesh.FaceRegn[ 2*iF+ilR ] = iR ;
    }
  }
  std::cout << "end-->build_EdgeFace" << std::endl << std::flush ;
}
void mesh2Dv_builder :: shrink_list( vector<int> & tmp_list ) {
  sort( tmp_list.begin(), tmp_list.end() ) ;
  int k = 0 ;
  for ( int i=1; i<tmp_list.size(); ++i ) {
    if ( tmp_list[i]!=tmp_list[k] ) {
      tmp_list[++k]=tmp_list[i] ;
    }
  }
  tmp_list.resize(k+1) ;
}
void mesh2Dv_builder :: build_boundary_lists() {
  std::cout << "begin build_boundary_lists" << std::endl << std::flush ;
  // introduce tmp vectors
  vector<int> tmp_vrtx ;
  vector<int> tmp_face ;
  vector<int> tmp_regn ;
  
  // gather all boundary items
  for ( int iF=0; iF<mesh.nF; ++iF ) {
    if ( mesh.FaceRegn[2*iF+1]==UNSET ) {
      tmp_face.push_back( iF ) ;
      tmp_regn.push_back( mesh.FaceRegn[2*iF+0] ) ;
      tmp_vrtx.push_back( mesh.FaceVrtx[nFV*iF+0] ) ;
      tmp_vrtx.push_back( mesh.FaceVrtx[nFV*iF+1] ) ;
    }
  }

  // shrink
  shrink_list( tmp_vrtx ) ;
  shrink_list( tmp_face ) ;
  shrink_list( tmp_regn ) ;

  // setup and copy
  mesh.bnd_vrtx.resize( tmp_vrtx.size() ) ;
  mesh.bnd_face.resize( tmp_face.size() ) ;
  mesh.bnd_regn.resize( tmp_regn.size() ) ;

  for ( int i=0; i<tmp_vrtx.size(); ++i ) { mesh.bnd_vrtx[i] = tmp_vrtx[i] ; }
  for ( int i=0; i<tmp_face.size(); ++i ) { mesh.bnd_face[i] = tmp_face[i] ; }
  for ( int i=0; i<tmp_regn.size(); ++i ) { mesh.bnd_regn[i] = tmp_regn[i] ; }

  std::cout << "end-->build_boundary_lists" << std::endl << std::flush ;
}

#include "build_face.hh"
#include "build_geom.hh"

// additional
#include "build_mesh.hh"

// DEBUG
#if 0
void mesh2Dv_builder :: print_vec_regn_face() {
  int nR = vec_regn_face.size() ;
  for ( int iR=0; iR<nR; ++iR ) {
    regn_face_struct & rfs = vec_regn_face[iR] ;
    int iF   = rfs.iGlb ;
    int ilF  = rfs.iloc ;
    int bval = rfs.bval ;
  } 
}
#endif
// fine DEBUG

#endif // end of -->> _MESH_2DV_BUILDER_HH
