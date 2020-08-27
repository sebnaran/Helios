#ifndef _MESH_2D_REFINER_HH
#define _MESH_2D_REFINER_HH

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// REGULAR SPLITTING of TRIANGLE-BASED MESH:
// - it splits an input triangle-based mesh by connecting edge midpoints of every triangle
// - the split mesh is guaranteed to satisfy the Delaunay condition only if the starting mesh 
//   is a Delaunay triangulation and its angles are (strictly) less that 90 degrees
// - external flags of regions  are propagated
// - external flags of vertices are propagated
// - new vertices from faces takes face flags
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void build_split_triangular_mesh( mesh_2Dv & mesh, mesh_2Dv & split_mesh ) {
  int nV = mesh.n_vertex()+mesh.n_face() ; // #vertices of split mesh
  int nR = 4*mesh.n_region() ;             // #regions  of split mesh
  vector<double> xV, yV ;
  vector<int>    regn_vlist ;
  vector<int>    fV, fR ;
  //-----------------------------------------------------------------------------------------
  // BUILD VERTEX DATA OF SPLIT MESH
  //-----------------------------------------------------------------------------------------
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    xV.push_back( mesh.coords_V(iV,0) ) ;
    yV.push_back( mesh.coords_V(iV,1) ) ;
    fV.push_back( mesh.get_fV(iV) ) ;
  }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    xV.push_back( mesh.coords_F(iF,0) ) ;
    yV.push_back( mesh.coords_F(iF,1) ) ;
    fV.push_back( mesh.get_fF(iF) ) ;
  }
  assert( nV==xV.size() && xV.size()==yV.size() ) ;
  //-----------------------------------------------------------------------------------------
  // BUILD REGION DATA OF SPLIT MESH
  // set region vertex lists
  int nv = mesh.n_vertex() ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> iv, ie ; // current lists of region vertices and edges
    mesh.get_regn_vrtx( iR, iv ) ;
    mesh.get_regn_face( iR, ie ) ;
    assert( iv.size() == 3 ) ;
    int regn_flag = mesh.get_fR( iR ) ;
    // first triangle
    regn_vlist.push_back(3) ;
    regn_vlist.push_back(   iv[0]) ;
    regn_vlist.push_back(nv+ie[0]) ;
    regn_vlist.push_back(nv+ie[2]) ;
    fR.push_back( regn_flag ) ;
    // second triangle
    regn_vlist.push_back(3) ;
    regn_vlist.push_back(   iv[1]) ;
    regn_vlist.push_back(nv+ie[1]) ;
    regn_vlist.push_back(nv+ie[0]) ;
    fR.push_back( regn_flag ) ;
    // third triangle
    regn_vlist.push_back(3) ;
    regn_vlist.push_back(   iv[2]) ;
    regn_vlist.push_back(nv+ie[2]) ;
    regn_vlist.push_back(nv+ie[1]) ;
    fR.push_back( regn_flag ) ;
    // fourth triangle
    regn_vlist.push_back(3) ;
    regn_vlist.push_back(nv+ie[0]) ;
    regn_vlist.push_back(nv+ie[1]) ;
    regn_vlist.push_back(nv+ie[2]) ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
  //-----------------------------------------------------------------------------------------
  // CALL BUILDER FOR SPLIT MESH (final construction)
  //-----------------------------------------------------------------------------------------
  mesh2Dv_builder mesh_builder( split_mesh ) ;
  mesh_builder . build_the_mesh( xV, yV, fV, regn_vlist, fR ) ; // !!!
}

#endif // end of _MESH_2D_REFINER_HH
