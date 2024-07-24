"""
    DirichletBoundaryCondition <: BoundaryCondition

Option to enforce that the surface plane `z=0` remains at temperature `T=0`
"""
struct DirichletBoundaryCondition <: BoundaryCondition end