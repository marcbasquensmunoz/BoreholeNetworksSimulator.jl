"""
    NeumannBoundaryCondition <: BoundaryCondition

Option to enforce zero heat flow at the surface plane `z=0`.
"""
struct NeumannBoundaryCondition <: BoundaryCondition end