using NamedDims

mutable struct feElementSet
    ChildElements
    Name::String
    feElementSet() = new()
end

mutable struct feSurfaceSet
    ChildElements
    ChildElements_LocalFace
    ChildNodes
    Name::String
    feSurfaceSet() = new()
end

mutable struct feNodeSet
    ChildNodes
    Variate
    Name::String
    feNodeSet() = new()
end

mutable struct feNode
    Basis::String
    ParentElements 
    Coordinates 
    isElementBoundaryNode::Bool
    isElementCornerNode::Bool
    isElementFaceNode::Bool
    isElementInternalNode::Bool
    isMeshBoundaryNode::Bool
    isMeshCornerNode::Bool
    isMeshFaceNode::Bool
    isMeshInternalNode::Bool
    ChildDOFS
    feNode() = new()
end

mutable struct feQuadrature
    Type::String
    Points
    Variate_Points
    Weights
    Variate_Weights
    Nₐ
    ∇Nₐ
    ∇ₓNₐ
    Jᵢⱼ
    feQuadrature() = new()
end

mutable struct feElement
    Basis::String
    BoundaryNodes::NamedDimsArray{(:local_node_id,)}
    ParentBlock::Int
    ChildNodes 
    CornerNodes 
    Degree 
    Dimension 
    NumVariates
    ElementFamily 
    FaceNodes 
    GlobalID 
    InternalNodes
    Quadrature::feQuadrature
    NumNodes::Int
    feElement() = new()
end

mutable struct MESH  
    Elements::NamedDimsArray{(:global_element_id,)}
    Nodes::NamedDimsArray{(:global_node_id,)}
    NodeSets::NamedDimsArray{(:global_nodeset_id,)}
    SurfaceSets::NamedDimsArray{(:global_surfaceset_id,)}
    MESH() = new()
end

struct ExodusElement
    ElementType::String
    ElementFaceOrder::NamedDimsArray{(:local_face_id,)}
    ElementNodeOrder::NamedDimsArray{(:local_node_id,)}
    FaceNodeOrder::NamedDimsArray{(:local_face_id,:local_node_id,)}
    isBoundaryNode::NamedDimsArray{(:local_node_id,)}
    isCornerNode::NamedDimsArray{(:local_node_id,)}
    isFaceNode::NamedDimsArray{(:local_node_id,)}
    isInternalNode::NamedDimsArray{(:local_node_id,)}
end