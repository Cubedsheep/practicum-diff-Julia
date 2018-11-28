include("solvers/f.jl")
include("solvers/Runge-Kutta.jl")
include("oef4.jl")

function testoef4(guess_y0, rounds, h_index)
    x0 = -11.998523280062
    z0 = 27.0
    guess_y0 = guess_y0
    guess_T = [3.0, 3.1]
    rounds = rounds

    gues1 = findy0(f, x0, z0, guess_y0, guess_T, rounds, h_index)
    #print(gues1[end])

    return gues1
end
