using LinearAlgebra


"""
berekent de afstand van de rechte opgespannen door de punten beginP en endP
tot het punt point
begin, end en point moeten 1d vectoren van dezelfde dimensie zijn
returns: Float64
"""
function dist_line_point(beginP::Array{Float64, 1}, endP::Array{Float64, 1},
     point::Array{Float64, 1});
    return norm(cross(endP-beginP, beginP-point))/ norm((endP-beginP))
end


"""
berekent de orthogonale projectie van het punt point op de rechte opgespannen door
beginP en endP
"""
function orth_projection(beginP, endP, point)
    P = point;
    P0= beginP;
    v = P-P0;
    s = endP - beginP;
    I = fill(1., 3, 3);
    projectie = (dot(v, s) / dot(s,s))*s;
    return projectie + P0
end

"""
berekent de afstand van het punt point tot het lijnstuk tussen corner1 en
corner 2. Returns de afstand tot het lijnstuk, en hoe ver de projectie
op het lijnstuk ligt als fractie van de lengte, tellend vanaf corner1.
"""
function dist_point_segment(corner1::Array{Float64, 1},
                            corner2::Array{Float64, 1},
                            point::Array{Float64, 1})
    # bereken de afstand tot de 2 hoekpunten en de afstand tot de rechte
    d_corner1 = norm(corner1 - point);
    d_corner2 = norm(corner2 - point);
    d_line = dist_line_point(corner1, corner2, point);

    # check of de orthogonale projectie op het lijnstuk ligt
    projection = orth_projection(corner1, corner2, point);
    # bereken de verhoudingen
    fractions = (corner1 - projection) ./ (corner1 - corner2);
    # test of alle verhoudingen tussen 0 en 1 liggen
    on_segment = all(x -> (x >= 0) && (x <= 1), fractions);
    fraction = fractions[1];
    if on_segment
        return (d_line, fraction);
    elseif fraction < 0
        return (d_corner1, 0.0);
    else
        return (d_corner2, 1);
    end
end

"""
berekent de volgende iteratiestap voor de methode van Runge-Kutta van orde 4
om de oplossing van de diff. vgl. X'=f(x) te benaderen met stapgrootte h
"""
function RK4_step(f::Function, Xi::Array{Float64, 1}, h::Float64)
    k1 = f(Xi);
    k2 = f(Xi + h*k1/2);
    k3 = f(Xi + h*k2/2);
    k4 = f(Xi + h*k3);
    return Xi + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

"""
berekent de eerste n iteratiestappen van de Runge-Kutta methode van orde 4
en returnt het tijdstip waarop, de kleinste afstand tot de beginwaarde bereikt
werd in het interval [n1, n2] en deze afstand
"""
function RK4dists(f::Function, X0::Array{Float64, 1}, h::Float64, n1::Int64,
                    n2::Int64)
    # voeg de eerste 2 punten toe aan de array met punten
    X1 = X0;
    X2 =  RK4_step(f, X1, h);

    # bereken de eerste n1 iteratiestappen
    for i = 2:n1
        X1 = X2;
        X2 = RK4_step(f, X1, h);
    end

    # bereken voor de volgende stappen ook steeds de afstand
    dist, t = dist_point_segment(X1, X2, X0)
    best = (dist, h*(n1+t-1))

    # itereer over de resterende stappen
    for i = (n1+1):n2
        X1 = X2;
        X2 = RK4_step(f, X1, h);
        dist, t = dist_point_segment(X1, X2, X0)
        if dist < best[1]
            best = (dist, h*(i+t-1))
        end
    end
    return best
end
