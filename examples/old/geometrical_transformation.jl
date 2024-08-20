function rotation(X,θ)
    rotation_matrix = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    return [GeometryTypes.Point2(rotation_matrix*x) for x in X]
end

function rotation_z(X,θ)
    rotation_matrix = [cos(θ) -sin(θ) 0.; sin(θ) cos(θ) 0.; 0. 0. 1]
    return [GeometryTypes.Point3(rotation_matrix*x) for x in X]
end