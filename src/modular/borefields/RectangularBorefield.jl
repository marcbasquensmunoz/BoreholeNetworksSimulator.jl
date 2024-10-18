
function RectangularBorefield(n1::Int, n2::Int, B1::T, B2::T, borehole_prototype::Borehole) where {T <: Number}
    positions = [((i-1)*B1, (j-1)*B2) for i in 1:n1 for j in 1:n2]
    EqualBoreholesBorefield(borehole_prototype = borehole_prototype, positions = positions)
end
