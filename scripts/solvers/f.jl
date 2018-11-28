function f(X::Array{Float64, 1})
    (x, y, z) = X;

    dx = -10.0*x+10.0*y;
    dy = 28.0*x-y-x*z;
    dz = -8/3*z+x*y;

    return [dx, dy, dz]
end
