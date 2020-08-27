// this is the 3-D implementation!
struct face_struct {
public:
  vector<int> vlist ; // original vertex list
  vector<int> slist ; // sorted   vertex list
  int iR  ;           // region of the face iF  
  int ilF ;           // local position of iF inside iR
  face_struct( vector<int> & _vlist, int _iR, int _ilF ) : 
    vlist(_vlist), slist(_vlist), iR(_iR), ilF(_ilF) {
    sort( slist.begin(), slist.end() ) ; 
  }
  ~face_struct() {}
} ;

bool operator<( const face_struct & F0, const face_struct & F1 ) {
  bool retval = F0.slist.size()<F1.slist.size() ;
  if ( F0.slist.size()==F1.slist.size() ) { 
    int i = 0 ;
    while( F0.slist[i]==F1.slist[i] && i<F0.slist.size() ) { ++i ; } 
    if ( i<F0.slist.size() ) { 
      retval = F0.slist[i]<F1.slist[i] ;
    } else {
      retval = F0.iR<F1.iR ;
    }
  }
  return retval ;
}

bool operator==( const face_struct & F0, const face_struct & F1 ) {
  bool retval(false) ;
  if ( F0.slist.size()==F1.slist.size() ) { 
    retval = true ;
    for ( int i=0; i<F0.slist.size() && retval; ++i ) {
      retval &= F0.slist[i]==F1.slist[i] ;
    }
  }
  return retval ;
}

bool operator!=( const face_struct & F0, const face_struct & F1 ) {
  return !( F0==F1 ) ;
}
// --- end of 3-D implementation

// INPUT:  regn_flist
// OUTPUT:
// (i)  build face_vlist from regn_flist
// (ii) collect face indices for mesh::RegnFace
void mesh2Dv_builder ::
build_face_vlist( vector< vector<regn_face_struct> > & vec_regn_face,   // output (construction)
		  vector<int>                        & face_vlist,      // output (new)
		  vector<int>                        & regn_vlist ) {   // input
  
  vector<int> vlist(2) ;
  vector<face_struct> vec_face_struct ;

  // ALL REGION FACES
  // fill vec_face_struct with instances of face_struct type
  // from all the region faces (any order is OK)
  int nR  = regn_vlist.back() ;
  int kR  = 0 ;                               // region pointer
  for ( int iR=0; iR<nR; ++iR ) {             // loop on all the regions
    int nRV = regn_vlist[kR++] ;              // #of vertices of region iR
    for ( int ilV=0; ilV<nRV; ++ilV ) {       // loop on the vertices of region iR
      vlist[0] = regn_vlist[kR+ilV        ] ; // first  vertex of the current face
      vlist[1] = regn_vlist[kR+(ilV+1)%nRV] ; // second vertex of the current face
      vec_face_struct.push_back( face_struct(vlist,iR,ilV) ) ;
    }
    kR += nRV ; // point to the next region in regn_vlist
  }

  // build an ORDERED set of faces represented by face_struct instances 
  // from all the region faces: internal faces are repetead twice
  sort( vec_face_struct.begin(), vec_face_struct.end() ) ;

  if ( false ) { // output for debugging
    int nF = vec_face_struct.size() ;
    for ( int k=0; k<nF; ++k ) { 
      std::cout << k << vec_face_struct[k].iR  << vec_face_struct[k].ilF ;
      std::cout << k << " -->> " ;
      int nFV = vec_face_struct[k].vlist.size() ; // #of vertices of face ilF
      for ( int ilV=0; ilV<nFV; ++ilV ) { // loop on the vertices of face ilF
	std::cout << " " << vec_face_struct[k].slist[ilV] ;
      }
      std::cout << "<<-- " ; 
      int iR  = vec_face_struct[k].iR  ; 
      int ilF = vec_face_struct[k].ilF ;
      std::cout << "\t" << " iR=" << iR << " ilF=" << ilF << " -->> " ;
      for ( int ilV=0; ilV<nFV; ++ilV ) { // loop on the vertices of face ilF
	std::cout << " " << vec_face_struct[k].vlist[ilV] ;
      }
      std::cout << std::endl ;
    }
  }

  { // make data structure face_vlist & set RegnFace
    // set the first face
    int iF = 0 ; // face index
    int kF = 0 ; // local pointer, be careful, kF==ilF is true only at this step!
    // --store the vertex list
    assert( vec_face_struct[kF].vlist.size()==2 ) ;
    for ( int ilV=0; ilV<vec_face_struct[kF].vlist.size(); ++ilV ) {
      face_vlist.push_back( vec_face_struct[kF].vlist[ilV] ) ;
    }
    { // --set RegnFace, face index is iF:=0
      int iR  = vec_face_struct[kF].iR  ;
      int ilF = vec_face_struct[kF].ilF ;
      //vec_regn_face[iR].push_back( regn_face_struct(iF,ilF,bool(true)) ) ;
      vec_regn_face[iR].push_back( regn_face_struct(iF,ilF) ) ;
    }
    // set the remaining faces
    for ( int il=1; il<vec_face_struct.size(); ++il ) {
      if ( vec_face_struct[il]!=vec_face_struct[kF] ) {
	// found a new face
	iF++    ; // update the face counter
	kF = il ; // reset  the face pointer
	// --store the vertex list
	assert( vec_face_struct[kF].vlist.size()==2 ) ;
	for ( int ilV=0; ilV<vec_face_struct[kF].vlist.size(); ++ilV ) {
	  face_vlist.push_back( vec_face_struct[kF].vlist[ilV] ) ;
	}
      }
      // --set RegnFace, face index is iF, ilF is the other instance of the face iF
      // --take the orientation of the face iF in the first region
      int  iR   = vec_face_struct[il].iR  ;
      int  ilF  = vec_face_struct[il].ilF ;
      vec_regn_face[iR].push_back( regn_face_struct(iF,ilF) ) ;
    }
    // finally, set the number of mesh faces
    face_vlist.push_back( iF+1 ) ;
  }

  if ( false ) {
    assert( nR == vec_regn_face.size() ) ;
    for ( int iR=0; iR<nR; ++iR ) {
      vector<regn_face_struct> & rfs_iR = vec_regn_face[iR] ;
      int nRF = rfs_iR.size() ;
      std::cout << "----------------------------------------" << std::endl ;
      std::cout << "iR=" << iR << " nRF=" << nRF << std::endl ; 
      for ( int ilF=0; ilF<nRF; ++ilF ) {
	std::cout << "rfs_iR[" << ilF << "].iGlb="  << rfs_iR[ilF].iGlb ;
	std::cout << " rfs_iR[" << ilF << "].iloc=" << rfs_iR[ilF].iloc ;
      }
    }
  } // end of -->> if ( false ) {...
}

// NEW NEW NEW
// INPUT:  
// face_vrtx_list
// face_regn_list
// OUTPUT:
// (i)  build face_vlist from regn_flist
// (ii) collect face indices for mesh::RegnFace
void mesh2Dv_builder ::
build_regn_flist( vector< vector<regn_face_struct> > & vec_regn_face,     // output (construction)
		  vector<int>                        & face_vrtx_list,    // input
		  vector<int>                        & face_regn_list ) { // input
    
  assert( face_vrtx_list.size()==face_regn_list.size() ) ;

  // get nR
  int nR = vec_regn_face.size() ;

  // build a list of faces for each region
  vector<int> regn_face_list[nR] ;
  {  //build regn_face_list ;
    int nF = face_regn_list.back() ;
    for ( int iF=0; iF<nF; ++iF ) {                 // loop on all the regions
      int iR0 = face_regn_list[2*iF+0] ; // first  region close to the current face
      int iR1 = face_regn_list[2*iF+1] ; // second region close to the current face
      if ( iR0!=UNSET ) { regn_face_list[iR0].push_back(iF) ; } 
      if ( iR1!=UNSET ) { regn_face_list[iR1].push_back(iF) ; } 
    }
  }
  
  if ( false ) { // printings for debug
    for ( int iR=0; iR<nR; ++iR ) {
      int nRF = regn_face_list[iR].size() ;
      std::cout << "----------------------------------------" << std::endl ;
      std::cout << "iR=" << iR << " nRF=" << nRF << std::endl ; 
      for ( int ilF=0; ilF<nRF; ++ilF ) {
	int iF = regn_face_list[iR][ilF] ;
	int iV0 = face_vrtx_list[ 2*iF   ] ;
	int iV1 = face_vrtx_list[ 2*iF+1 ] ;
	std::cout << "ilF=" << ilF << " iF=" << iF << " iV0=" << iV0 << " iV1=" << iV1 << std::endl ;
      }
    }
    exit(0) ;
  }

  int nF = face_regn_list.back() ;
  vector<bool> face_set(nF,true) ;

  for ( int iR=0; iR<nR; ++iR ) {

    int nRF = regn_face_list[iR].size() ;
    vector<bool> fmask(nRF,true) ;
    vector<int>  flist(nRF) ;
    vector<int>  fbool(nRF) ;

    // insert the starting face
    int iF  = regn_face_list[iR][0] ;
    int iV0 = face_vrtx_list[ 2*iF   ] ;
    int iV1 = face_vrtx_list[ 2*iF+1 ] ;

    fmask[0] = false ;
    flist[0] = iF ;

    int last_iV(UNSET) ;
    if ( iR==face_regn_list[2*iF] ) {
      last_iV = iV1 ;
      fbool[0] = true ;
    } else {
      last_iV = iV0 ;
      fbool[0] = false ;
    }
    assert( last_iV!=UNSET ) ;

    int ikF = 1 ; // the first face has already been inserted
    while ( ikF<nRF ) { // repeat the loop until there are faces to be processed

      for ( int ilF=1; ilF<nRF; ++ilF ) {	
	if ( fmask[ilF] ) {                   // face not yet found
	  int iF  = regn_face_list[iR][ilF]  ;
	  int iV0 = face_vrtx_list[ 2*iF   ] ;
	  int iV1 = face_vrtx_list[ 2*iF+1 ] ;
	  //std::cout << "searching --> " << " iF=" << iF << " iV0=" << iV0 << " iV1=" << iV1 << std::endl ;
	  if ( iV0==last_iV ) {         // found iV0
	    fmask[ilF] = false ;
	    flist[ikF] = iF ;
	    fbool[ikF] = true ;
	    last_iV    = iV1 ;
	    ++ikF ;
	    //std::cout << "found -->>" << " iF=" << iF << " last_iV=" << last_iV << std::endl ;
	    break ;
	  } else if ( iV1==last_iV ) {  // found iV1
	    fmask[ilF] = false ;
	    flist[ikF] = iF ;
	    fbool[ikF] = false ;
	    last_iV    = iV0 ;
	    ++ikF ;
	    //std::cout << "found -->>" << " iF=" << iF << " last_iV=" << last_iV << std::endl ;
	    break ;
	  } // end of -->> if ( iV0==last_iV ) {...
	} // end of -->> if ( fmask[ilF] ) {...
      } // end of -->> for ( int ilF=0; ilF<nRF; ++ilF ) {...
    } // end of -->> while ( ikF<nRF ) {...

    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF   = flist[ilF] ;
      vec_regn_face[iR].push_back( regn_face_struct(iF,ilF) ) ;
    }
  }

  assert( nR == vec_regn_face.size() ) ;
  if ( false ) {
    for ( int iR=0; iR<nR; ++iR ) {
      vector<regn_face_struct> & rfs_iR = vec_regn_face[iR] ;
      int nRF = rfs_iR.size() ;
      std::cout << "----------------------------------------" << std::endl ;
      std::cout << "iR=" << iR << " nRF=" << nRF << std::endl ; 
      for ( int ilF=0; ilF<nRF; ++ilF ) {
	std::cout << "rfs_iR[" << ilF << "].iGlb="  << rfs_iR[ilF].iGlb ;
	std::cout << " rfs_iR[" << ilF << "].iloc=" << rfs_iR[ilF].iloc ;
      }
    }
  }
}

void mesh2Dv_builder::build_the_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, 
				      vector<int> & face_vrtx_list, vector<int> & face_regn_list, 
				      vector<int> & fF ) {
  std::cout << "begin  ---> mesh2Dv_builder::build_the_mesh using face datasets <----" << std::endl << std::flush ;

  // check consistency
  assert( face_vrtx_list.size()==face_regn_list.size() ) ;
  assert( face_vrtx_list.size()%2==1 ) ;

  // -- set coordinates in mesh
  set_coords( xV, yV ) ;

  // determine nR 
  int nR = UNSET ;
  { 
    int nF = face_vrtx_list.back() ;
    for ( int iF=0; iF<nF; ++iF ) {    // loop on all the regions
      int iR0 = face_regn_list[2*iF+0] ; // first  region close to the current face
      int iR1 = face_regn_list[2*iF+1] ; // second region close to the current face
      nR = max( nR, max( iR0, iR1 ) ) ;
    }
  }
  nR += 1 ; 
  
  // -- core for building primary dataset
  vector< vector<regn_face_struct> > vec_regn_face(nR) ;
  vector<int> face_vlist ;
  build_regn_flist( vec_regn_face, face_vrtx_list, face_regn_list ) ;
  build_RegnFace( vec_regn_face, nR ) ;
  
  // set face-vertex data structure
  build_FaceVrtx( face_vrtx_list ) ;
  
  // -- set external flags for regions, faces and vertices
  mesh.fR.resize(nR) ;
  for ( int iR=0; iR<nR; ++iR ) { mesh.fR[iR] = UNSET ; }
  // ---
  int nF = face_vrtx_list.back() ;
  mesh.fF.resize(nF) ;
  for ( int iF=0; iF<nF; ++iF ) { mesh.fF[iF] = fF[iF] ; }
  // ---
  int nV = fV.size() ;
  mesh.fV.resize(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { mesh.fV[iV] = fV[iV] ; }  

  // -- tranpose datasets
  build_VrtxFace() ;
  build_FaceRegn() ;

  // -- build boundary lists
  build_boundary_lists() ;
  
  // -- compute/set last geometric quantities
  setup_geom_factors() ;
  
  // final reset of logical flags
  mesh.reset_status_flags() ;
  std::cout << "end of ---> mesh2Dv_builder::build_the_mesh using face datasets <----" << std::endl << std::flush ;
}

void mesh2Dv_builder::reorder_faces( vector<int> & face_perm ) {

  // vectors
  vector<int> FaceVrtx ;
  vector<int> FaceRegn ;

  // lists of boundary items
  vector<int> bnd_face ;

  // external flags
  vector<int> fF ;

  // mesh size
  int nV = mesh.n_vertex() ;
  int nF = mesh.n_face  () ;
  int nR = mesh.n_region() ;

  // permute RegnFace
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = mesh.n_regn_face(iR) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      mesh.RegnFace(iR,ilF) = face_perm[ mesh.regn_face(iR,ilF) ] ;
    }
  }

  // permute VrtxFace
  for ( int iV=0; iV<nV; ++iV ) {
    int nVF = mesh.n_vrtx_face(iV) ;
    for ( int ilF=0; ilF<nVF; ++ilF ) {
      mesh.VrtxFace(iV,ilF) = face_perm[ mesh.vrtx_face(iV,ilF) ] ;
    }
  }

  // permute FaceVrtx
  FaceVrtx.resize( 2*nF ) ;
  assert( FaceVrtx.size()==mesh.FaceVrtx.size() ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int new_iF = face_perm[ iF ] ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      FaceVrtx[nFV*new_iF+ilV] = mesh.face_vrtx(iF,ilV) ;
    }
  }
  for ( int i=0; i<FaceVrtx.size(); ++i ) {
    mesh.FaceVrtx[i] = FaceVrtx[i] ;
  }

  // permute FaceRegn ;
  int nFR = 2 ;
  FaceRegn.resize( nFR*nF ) ;
  assert( FaceRegn.size()==mesh.FaceRegn.size() ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int new_iF = face_perm[ iF ] ;
    for ( int ilR=0; ilR<nFR; ++ilR ) {
      FaceRegn[nFR*new_iF+ilR] = mesh.face_regn(iF,ilR) ;
    }
  }
  for ( int i=0; i<FaceRegn.size(); ++i ) {
    mesh.FaceRegn[i] = FaceRegn[i] ;
  }

  // permute bnd_face ;
  int nbF = mesh.n_bface() ;
  bnd_face.resize( nbF ) ;
  for ( int ilF=0; ilF<nbF; ++ilF ) {
    bnd_face[ilF] = face_perm[ mesh.get_bnd_face(ilF) ] ;
  }
  for ( int ilF=0; ilF<nbF; ++ilF ) {
    mesh.bnd_face[ilF] = bnd_face[ilF] ;
  }

  // permute fF
  fF.resize( nF ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int new_iF = face_perm[ iF ] ;
    fF[new_iF] = mesh.get_fF( iF ) ;
  }
  for ( int iF=0; iF<nF; ++iF ) {
    mesh.fF[iF] = fF[iF] ;
  }

  // eventually...
  recompute_geom_factors() ;
}
