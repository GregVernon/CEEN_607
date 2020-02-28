using NamedDims

mutable struct feElementSet
    ChildElements
    Dimension 
    NumVariates
    ElementFamily 
    LC_Type
    LC_Direction
    LC_Magnitude
    BC_Type
    BC_DOF
    BC_Value
    Name::String
    feElementSet() = new()
end

mutable struct feSurfaceSet
    ChildElements
    ChildElements_LocalFace
    ChildNodes
    LC_Type
    LC_Direction
    LC_Magnitude
    BC_Type
    BC_DOF
    BC_Value
    Name::String
    feSurfaceSet() = new()
end

mutable struct feNodeSet
    ChildNodes
    LC_Type
    LC_Direction
    LC_Magnitude
    BC_Type
    BC_DOF
    BC_Value
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

mutable struct feQuadraturePoint
    feQuadraturePoint() = new()
    Coordinates
    Weights
    Nₐ
    ∇Nₐ
    ∇ₓNₐ
    Jᵢⱼ
end

mutable struct feQuadrature
    Type::String
    QuadraturePoints::NamedDimsArray{(:local_qp_id,)}
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
    SideNodes
    GlobalID 
    InternalNodes
    Quadrature::NamedDimsArray{(:local_side_id,)}
    NumNodes::Int
    feElement() = new()
end

mutable struct MESH  
    Elements::NamedDimsArray{(:global_element_id,)}
    Nodes::NamedDimsArray{(:global_node_id,)}
    NodeSets::NamedDimsArray{(:global_nodeset_id,)}
    SurfaceSets::NamedDimsArray{(:global_surfaceset_id,)}
    ElementSets::NamedDimsArray{(:global_elementset_id,)}
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

################### Input Parameters ###################
mutable struct InputParams
    BoundaryConditions
    LoadConditions
    InputParams() = new()
end

mutable struct GeometryCard
    Filename
    GeometryCard() = new()
end

mutable struct BoundaryCondition
    Type
    NodeSetName
    DOF
    Value
    BoundaryCondition() = new()
end

mutable struct LoadCondition
    Type
    ElementSetName
    SurfaceSetName
    Direction
    Magnitude
    LoadCondition() = new()
end

mutable struct BodyCondition
    Type
    ElementSetName
    Value
    BodyCondition() = new()
end

