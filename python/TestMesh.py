from MeshHelios import HeliosMesh

#This is a simple test to check that we can cosntruct HeliosMeshes
def test_MeshHeliosInit():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    BoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations)
    assert (Nodes == TestMesh.Nodes and EdgeNodes == TestMesh.EdgeNodes and ElementEdges == TestMesh.ElementEdges)

def test_internalNodes():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    BoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations)
    #TestMesh.MakeIntNodes()
    assert (TestMesh.InternalNodes == [4])

def test_MakeDictionaries():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    BoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]

    TestMesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations)
    
    EdgestoCells = [[1],[0,1],[3],[2,3],[3],[2],[2],[1],[0],[0],[1,2],[0,3]]
    NodestoCells = [[1],[0,1],[0],[1,2],[0,1,2,3],[0,3],[2],[2,3],[3]]
    #TestMesh.MakeDictionaries()
    assert (TestMesh.EdgestoCells == EdgestoCells and TestMesh.NodestoCells == NodestoCells)