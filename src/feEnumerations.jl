module feEnumerations

export enum_BoundaryCondition
export enum_LoadCondition
export enum_MeshEntity
export enum_BasisFunction

@enum enum_BoundaryCondition dirichlet neumann
@enum enum_LoadCondition pressure traction force
@enum enum_MeshEntity element node
@enum enum_BasisFunction bernstein lagrange legendre 

end