
get_Nb(X) = Int(size(X)[1] / 4)
extract_Tfin(X) = @view X[1:2:2*get_Nb(X), :]
extract_Tfout(X) = @view X[2:2:2*get_Nb(X), :]
extract_Tb(X) = @view X[2*get_Nb(X)+1:3*get_Nb(X), :]
extract_q(X) = @view X[3*get_Nb(X)+1:end, :]