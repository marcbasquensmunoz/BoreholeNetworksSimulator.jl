
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
#setup(::MidPointApproximation, borefield::Borefield, i, j) = get_stp(borefield, i, j)
