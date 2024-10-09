
image(s::SegmentToSegment) = SegmentToSegment(D1=-s.D1, H1=-s.H1, D2=s.D2, H2=s.H2, σ=s.σ)
image(s::SegmentToPoint) = SegmentToPoint(D=-s.D-s.H, H=s.H, z=s.z, σ=s.σ)

function get_sts(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    SegmentToSegment(D1=Di, H1=Hi, D2=Dj, H2=Hj, σ=σ)
end

function get_stp(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    SegmentToPoint(σ = σ, D = Di, H = Hi, z = Dj + Hj/2)
end

setup(::MeanApproximation, borefield::Borefield, i, j) = get_sts(borefield, i, j)
setup(::MidPointApproximation, borefield::Borefield, i, j) = get_stp(borefield, i, j)

setup_type(::MeanApproximation) = SegmentToSegment{Float64}
setup_type(::MidPointApproximation) = SegmentToPoint{Float64}