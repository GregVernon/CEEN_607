using NamedDims

mutable struct feElementSet
    ChildElements::Array{Int,1}
    Dimension::Int
    NumVariates::Int
    ElementFamily::String
    LC_Type::NamedDimsArray{(:lc_id,), feEnumerate.enum_LoadCondition}           # push(NS.LC_Type, feEnumerate.pressure)
    LC_Direction::NamedDimsArray{(:lc_id,),Array{Float64,1}}                     # NS.LC_Direction = [NS.LC_Direction; [1. 2.] ]
    LC_Magnitude::NamedDimsArray{(:lc_id,), Float64}                             # push!(NS.LC_Magnitude, 10.) or NS.LC_Magnitude = [NS.LC_Magnitude; [10.]]
    LC_NL_Solve_Scale::Float64
    BC_Type::NamedDimsArray{(:bc_id,), feEnumerate.enum_BoundaryCondition}
    BC_DOF::NamedDimsArray{(:bc_id,), Int}
    BC_Value::NamedDimsArray{(:bc_id,), Float64}
    BC_NL_Solve_Scale::Float64
    Name::String
    feElementSet() = new()
end

mutable struct feSurfaceSet
    ChildElements::Array{Int,1}
    ChildElements_LocalFace::Array{Int,1}
    ChildNodes::Array{Int,1}
    LC_Type::NamedDimsArray{(:lc_id,), feEnumerate.enum_LoadCondition}           # push(NS.LC_Type, feEnumerate.pressure)
    LC_Direction::NamedDimsArray{(:lc_id,),Array{Float64,1}}                     # NS.LC_Direction = [NS.LC_Direction; [1. 2.] ]
    LC_Magnitude::NamedDimsArray{(:lc_id,), Float64}                             # push!(NS.LC_Magnitude, 10.) or NS.LC_Magnitude = [NS.LC_Magnitude; [10.]]
    LC_NL_Solve_Scale::Float64
    BC_Type::NamedDimsArray{(:bc_id,), feEnumerate.enum_BoundaryCondition}
    BC_DOF::NamedDimsArray{(:bc_id,), Int}
    BC_Value::NamedDimsArray{(:bc_id,), Float64}
    BC_NL_Solve_Scale::Float64
    Name::String
    feSurfaceSet() = new()
end

mutable struct feNodeSet
    ChildNodes::Array{Int,1}
    LC_Type::NamedDimsArray{(:lc_id,), feEnumerate.enum_BoundaryCondition}           # push!(NS.LC_Type, feEnumerate.pressure)
    LC_Direction::NamedDimsArray{(:lc_id,),Array{Float64,1}}                     # NS.LC_Direction = [NS.LC_Direction; [1. 2.] ]
    LC_Magnitude::NamedDimsArray{(:lc_id,), Float64}                             # push!(NS.LC_Magnitude, 10.) or NS.LC_Magnitude = [NS.LC_Magnitude; [10.]]
    LC_NL_Solve_Scale::Float64
    BC_Type::NamedDimsArray{(:bc_id,), feEnumerate.enum_BoundaryCondition}
    BC_DOF::NamedDimsArray{(:bc_id,), Int}
    BC_Value::NamedDimsArray{(:bc_id,), Float64}
    BC_NL_Solve_Scale::Float64
    Name::String
    feNodeSet() = new()
end

mutable struct feNode
    Basis::String
    ParentElements::Array{Int,1}
    Coordinates::NamedDimsArray{(:ℝᴺ,)}
    ChildDOFS::NamedDimsArray{(:local_dof_id,)}
    feNode() = new()
end

mutable struct feQuadraturePoint
    ℙ::NamedDimsArray{(:ℙᴺ,)}  # 1xN vector
    𝒲::Float64
    𝓝::NamedDimsArray{(:local_node_id,)}
    ∇𝓝::NamedDimsArray{(:local_node_id, :ℙᴺ,)}
    ∇ₓ𝓝::NamedDimsArray{(:local_node_id, :ℝᴺ)}
    Jᵢⱼ::NamedDimsArray{(:ℝᴺ,:ℙᴺ)}
    α::Float64
    ñ::NamedDimsArray{(:ℙᴺ,)}  # 1xN vector
    feQuadraturePoint() = new()
end

mutable struct feQuadrature
    Type::String
    QuadraturePoints::NamedDimsArray{(:local_qp_id,)}
    feQuadrature() = new()
end

mutable struct feElement
    Basis::String
    BoundaryNodes::NamedDimsArray{(:local_node_id,)}
    ChildNodes::Array{Int,1}
    Degree::Array{Int,1}
    ElementFamily::String
    GlobalID::Int
    ParentBlock::Int
    SideNodes::NamedDimsArray{(:local_side_id,)}
    StiffnessMatrix::NamedDimsArray{(:local_dof_id, :local_dof_id,)}
    ExternalForceVector::NamedDimsArray{(:local_dof_id,)}
    InternalForceVector::NamedDimsArray{(:local_dof_id,)}
    Quadrature::NamedDimsArray{(:local_side_id,)}
    NumNodes::Int
    𝓝
    ∇𝓝
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

