module feEnumerate
export enum_BoundaryCondition
export enum_LoadCondition

@enum enum_BoundaryCondition unconstrained dirichlet neumann
@enum enum_LoadCondition noload body pressure traction force

end