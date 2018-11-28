include("solvers/Runge-Kutta.jl")
using BenchmarkTools


"""
het rechterlid van de diff vgl X' = f(X)
"""
function f(X::Array{Float64, 1})
    (x, y, z) = X;

    dx = -10.0*x+10.0*y;
    dy = 28.0*x-y-x*z;
    dz = -8/3*z+x*y;

    return [dx, dy, dz]
end

X0 = [-13.763610682134, -19.578751942452, 27.0];
h = 1e-4;
T1 = 1.55; T2 = 1.56;
n1 = floor(Int64, T1/h);
n2 = ceil(Int64, T2/h);

closest = RK4dists(f, X0, h, n1, n2);
println(closest)

#println(@btime RK4_dists($f, $X0, $h, $n1, $n2) )
