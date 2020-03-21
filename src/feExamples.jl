using NamedDims
include("feDatastructure.jl")

############################################
############## 1 x 1 Element ###############
############################################

ELEM = [feElement()]
ELEM[1].Basis = "Lagrange"
ELEM[1].BoundaryNodes = NamedDimsArray{(:local_node_id,)}([1,2,3,4])
ELEM[1].ParentBlock = 1
ELEM[1].ChildNodes = NamedDimsArray{(:global_node_id,)}([1,2,3,4])
ELEM[1].CornerNodes = NamedDimsArray{(:local_node_id,)}([1,2,3,4])
ELEM[1].Degree = [1,1]
ELEM[1].Dimension = 2
ELEM[1].NumVariates = 2
ELEM[1].ElementFamily = "QUAD"
ELEM[1].FaceNodes = NamedDimsArray{(:local_face_id,)}([NamedDimsArray{(:local_node_id,)}([1,2,3,4])])
ELEM[1].GlobalID = 1
ELEM[1].InternalNodes = NamedDimsArray{(:local_node_id,)}([])
ELEM[1].NumNodes = 4

ELEM[1].Quadrature = feQuadrature()
ELEM[1].Quadrature.Type = "Gauss-Legendre"
ELEM[1].Quadrature.Points = NamedDimsArray{(:local_quadrature_id,:ℝᴺ,)}([-1/√3 -1/√3; +1/√3 -1/√3;-1/√3 +1/√3; +1/√3 +1/√3])
ELEM[1].Quadrature.Variate_Points = NamedDimsArray{(:global_variate_id,)}([NamedDimsArray{(:variate_quadrature_id,)}([-1/√3, +1/√3]), NamedDimsArray{(:variate_quadrature_id,)}([-1/√3, +1/√3])])
ELEM[1].Quadrature.Weights = NamedDimsArray{(:local_quadrature_id,)}([1, 1, 1, 1])
ELEM[1].Quadrature.Variate_Weights = NamedDimsArray{(:global_variate_id,)}([NamedDimsArray{(:variate_quadrature_id,)}([1, 1]), NamedDimsArray{(:variate_quadrature_id,)}([1, 1])])
ELEM[1].Quadrature.Nₐ = NamedDimsArray{(:local_quadrature_id,)}([((1-(-1/√3))/2) * ((1-(-1/√3))/2), ((1-(-1/√3))/2) * ((1+(-1/√3))/2), ((1-(-1/√3))/2) * ((1-(-1/√3))/2), ((1+(-1/√3))/2) * ((1+(-1/√3))/2)])
ELEM[1].Quadrature.∇Nₐ = NamedDimsArray{(:local_quadrature_id,:ℙᴺ,)}([-1/2 -1/2; 1/2 -1/2; -1/2 1/2; 1/2 1/2])
ELEM[1].Quadrature.∇ₓNₐ = NamedDimsArray{(:local_quadrature_id,:ℝᴺ,)}([-1/2 -1/2; 1/2 -1/2; -1/2 1/2; 1/2 1/2])
ELEM[1].Quadrature.Jᵢⱼ = NamedDimsArray{(:ℝᴺ, :ℙᴺ,)}([1.0 0.0; 0.0 1.0])


NODE = [feNode(), feNode(), feNode(), feNode()]
NODE[1].Basis = "Lagrange"
NODE[1].ParentElements = NamedDimsArray{(:global_elem_id,)}([1])
NODE[1].Coordinates = NamedDimsArray{(:ℝᴺ,)}([-1.0, -1.0])
NODE[1].isElementBoundaryNode = true
NODE[1].isElementCornerNode = true
NODE[1].isElementFaceNode = true
NODE[1].isElementInternalNode = false
NODE[1].isMeshBoundaryNode = true
NODE[1].isMeshCornerNode = true
NODE[1].isMeshFaceNode = true
NODE[1].isMeshInternalNode = false
NODE[1].ChildDOFS = NamedDimsArray{(:global_dof_id,)}([1,2])

NODE[2].Basis = "Lagrange"
NODE[2].ParentElements = NamedDimsArray{(:global_elem_id,)}([1])
NODE[2].Coordinates = NamedDimsArray{(:ℝᴺ,)}([+1.0, -1.0])
NODE[2].isElementBoundaryNode = true
NODE[2].isElementCornerNode = true
NODE[2].isElementFaceNode = true
NODE[2].isElementInternalNode = false
NODE[2].isMeshBoundaryNode = true
NODE[2].isMeshCornerNode = true
NODE[2].isMeshFaceNode = true
NODE[2].isMeshInternalNode = false
NODE[2].ChildDOFS = NamedDimsArray{(:global_dof_id,)}([3,4])

NODE[3].Basis = "Lagrange"
NODE[3].ParentElements = NamedDimsArray{(:global_elem_id,)}([1])
NODE[3].Coordinates = NamedDimsArray{(:ℝᴺ,)}([-1.0, +1.0])
NODE[3].isElementBoundaryNode = true
NODE[3].isElementCornerNode = true
NODE[3].isElementFaceNode = true
NODE[3].isElementInternalNode = false
NODE[3].isMeshBoundaryNode = true
NODE[3].isMeshCornerNode = true
NODE[3].isMeshFaceNode = true
NODE[3].isMeshInternalNode = false
NODE[3].ChildDOFS = NamedDimsArray{(:global_dof_id,)}([5,6])

NODE[4].Basis = "Lagrange"
NODE[4].ParentElements = NamedDimsArray{(:global_elem_id,)}([1])
NODE[4].Coordinates = NamedDimsArray{(:ℝᴺ,)}([+1.0, +1.0])
NODE[4].isElementBoundaryNode = true
NODE[4].isElementCornerNode = true
NODE[4].isElementFaceNode = true
NODE[4].isElementInternalNode = false
NODE[4].isMeshBoundaryNode = true
NODE[4].isMeshCornerNode = true
NODE[4].isMeshFaceNode = true
NODE[4].isMeshInternalNode = false
NODE[4].ChildDOFS = NamedDimsArray{(:global_dof_id,)}([7,8])

SET_ELEM = [feElementSet()]
SET_ELEM[1].ChildElements=NamedDimsArray{(:global_elem_id,)}([1])
SET_ELEM[1].Name="Domain"

SET_SURF = [feSurfaceSet()]
SET_SURF[1].ChildElements = NamedDimsArray{(:global_elem_id,)}([1])
SET_SURF[1].ChildElements_LocalFace = NamedDimsArray{(:local_face_id,)}([2])
SET_SURF[1].ChildNodes = NamedDimsArray{(:global_node_id,)}([2, 4])
SET_SURF[1].Name = "Load Surface"

SET_NODE = [feNodeSet()]
SET_NODE[1].ChildNodes = NamedDimsArray{(:global_node_id,)}([1,3])
SET_NODE[1].Variate = NamedDimsArray{(:global_variate_id,)}([1,2])
SET_NODE[1].Name = "Hold Nodes"

# Global arrays
eCONN = NamedDimsArray{(:local_node_id,:global_elem_id,)}(collect(transpose([1 2 3 4])))
nodeCoord = NamedDimsArray{(:local_node_id,:ℝᴺ,)}([-1.0 -1.0; +1.0 -1.0; -1.0 +1.0; +1.0 +1.0])


############################################
############## 2 x 2 Element ###############
############################################

