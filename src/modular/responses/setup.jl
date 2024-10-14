
image(s::SegmentToSegment) = SegmentToSegment(D1=-s.D1, H1=-s.H1, D2=s.D2, H2=s.H2, σ=s.σ)
image(s::SegmentToPoint) = SegmentToPoint(D=-s.D-s.H, H=s.H, z=s.z, σ=s.σ)
image(s::MovingSegmentToSegment) = MovingSegmentToSegment(x=s.x, y=s.y,  v=s.v, D1=-s.D1, H1=-s.H1, D2=s.D2, H2=s.H2)
image(s::MovingSegmentToPoint) = MovingSegmentToPoint(x=s.x, y=s.y, z=s.z, D=-s.D-s.H, H=s.H, v=s.v)

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

function get_msts(borefield::Borefield, i, j, v)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    x = i == j ? 0. : xj - xi
    y = i == j ? get_rb(borefield, i) : yj - yi
    MovingSegmentToSegment(x = x, y = y, v = v, D1=Di, H1=Hi, D2=Dj, H2=Hj)
end

function get_mstp(borefield::Borefield, i, j, v)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    x = i == j ? 0. : xj - xi
    y = i == j ? get_rb(borefield, i) : yj - yi
    MovingSegmentToPoint(x = x, y = y, z = Dj + Hj/2, v = v, D = Di, H = Hi)
end


setup(::MeanApproximation, ::GroundMedium, borefield::Borefield, i, j) = get_sts(borefield, i, j)
setup(::MidPointApproximation, ::GroundMedium, borefield::Borefield, i, j) = get_stp(borefield, i, j)

setup(::MeanApproximation, medium::FlowInPorousMedium, borefield::Borefield, i, j) = get_msts(borefield, i, j, medium.vt)
setup(::MidPointApproximation, medium::FlowInPorousMedium, borefield::Borefield, i, j) = get_mstp(borefield, i, j, medium.vt)

setup_type(::MeanApproximation, ::GroundMedium) = SegmentToSegment{Float64}
setup_type(::MidPointApproximation, ::GroundMedium) = SegmentToPoint{Float64}
setup_type(::MeanApproximation, ::FlowInPorousMedium) = MovingSegmentToSegment{Float64}
setup_type(::MidPointApproximation, ::FlowInPorousMedium) = MovingSegmentToPoint{Float64}