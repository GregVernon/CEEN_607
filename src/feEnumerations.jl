module feEnumerations

export enum_BoundaryCondition
export enum_LoadCondition

@enum enum_BoundaryCondition dirichlet neumann
@enum enum_LoadCondition pressure traction force
@enum enum_MeshEntity element node

end