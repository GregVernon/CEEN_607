module feDatastruct

export feElement, feNode

mutable struct feMesh
    Elements # DONE
    Nodes # DONE
    ElementSets
    SurfaceSets
    NodeSets
end


mutable struct feElement
    Basis
    BoundaryNodes # DONE
    ParentBlocks # DONE
    ChildNodes # DONE
    CornerNodes # DONE
    Degree # DONE
    Dimension # DONE
    ElementFamily # DONE
    FaceNodes # DONE
    GlobalID # DONE
    InternalNodes # DONE
    Quadrature
    NumNodes # DONE
    feElement() = new()
end

mutable struct feNode
    Basis
    ParentElements # DONE
    Coordinates # DONE
    isElementBoundaryNode # DONE
    isElementCornerNode # DONE
    isElementFaceNode # DONE
    isElementInternalNode # DONE
    isMeshBoundaryNode
    isMeshCornerNode
    isMeshFaceNode
    isMeshInternalNode
    ChildDOFS # DONE
    feNode() = new()
end

mutable struct feElementSet
    ChildElements
end

mutable struct feSurfaceSet
    Child
end

mutable struct feNodeSet
    ChildNodes
end

mutable struct Quadrature
    Type
    Points
    Weights
    Quadrature() = new()
end

struct ExodusElement
    ElementType # DONE
    ElementNodeOrder # DONE
    FaceNodeOrder # DONE
    EdgeNodeOrder # DONE
    isBoundaryNode # DONE
    isCornerNode # DONE
    isFaceNode # DONE
    isInternalNode # DONE
end

function makeExodusElement(elem_type)
    ElementType = elem_type
    if elem_type == "BAR2"
        ElementNodeOrder = [1 2]
        FaceNodeOrder = [[1 2]]
        isBoundaryNode = fill(true,2)
        isCornerNode = fill(true,2)
        isFaceNode = fill(true,2)
        isInternalNode = fill(false,2)
    end

    if elem_type == "QUAD4"
        ElementNodeOrder = [1 2 4 3]
        FaceNodeOrder = [[1 2 4 3]]
        isBoundaryNode = fill(true,4)
        isCornerNode = fill(true,4)
        isFaceNode = fill(true,4)
        isInternalNode = fill(false,4)
    end

    if elem_type == "HEX8"
        ElementNodeOrder = [1 2 4 3 5 6 8 7]
        ElementFaceOrder = [4 2 1 3 5 6]
        FaceNodeOrder = [[1 2 6 5],[2 4 8 6],[4 3 7 8],[1 5 7 3],[1 3 4 2],[5 6 8 7],[],[]]
        isBoundaryNode = fill(true,8)
        isCornerNode = fill(true,8)
        isFaceNode = fill(true,8)
        isInternalNode = fill(false,8)
    end
    GE = ExodusElement(ElementType,ElementNodeOrder,ElementFaceOrder,FaceNodeOrder,isBoundaryNode,isCornerNode,isFaceNode,isInternalNode)
    return GE
end


end