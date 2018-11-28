@everywhere begin
    include("solvers/Runge-Kutta.jl")
    include("solvers/f.jl")

end

"""
initialiseerd de arrays die aanduiden in welk gebied de kleinste
afstand gezocht wordt, in termen van de index
"""
function init_N(h, guess_T, num)
    H = fill(h, num);
    n1 = floor(Int64, guess_T[1]/h); n2 = ceil(Int64, guess_T[2]/h);
    N1 = fill(n1, num); N2 = fill(n2, num);
    return (H, n1, n2, N1, N2)
end

"""
functie om beginvoorwaarden te zoeken die een periodieke oplossing geven.
Zoekt in het interval guess_y0 naar een waarde y0 zodat (x0, y0, z0) een
beginvoorwaarde van een periodieke oplossing is met periode T in guess_T.
Het verfijnt deze benadering rounds keer.
"""
function findy0(f, x0, z0, guess_y0, guess_T, rounds, h_index)
    # arrays die info bijhouden over de haalbare precisie met
    # een bepaalde stapgrootte, voor verkleinen stapgrootte indien
    # nodig.
    num_h = 4;
    delta_arr = [1e-5, 1e-7, 1e-9, 1e-11];
    h_arr = [1e-4, 1e-5, 1e-6, 1e-7];
    h = h_arr[h_index];

    # initializeren variabelen
    # vul arrays met de stapgrootte, aantal stappen en rechterlid
    # om parallel de vergelijkingen
    # op te lossen
    H, n1, n2, N1, N2 = init_N(h, guess_T, 7);
    F = fill(f, 7);

    # los de vergelijkingen op en zoek de beginvoorwaarde die de beste
    # benadering geeft
    X = [[x0, guess_y0[1] + i*(guess_y0[2]-guess_y0[1])/6, z0] for i=0:6];
    dists = pmap(RK4dists, F, X, H, N1, N2);
    # array met afstanden
    d = [dist[1] for dist in dists];
    # zoek de index van de kleinste afstand
    mini = minimum(d);
    index = findfirst(isequal(mini), d);

    # als de beste waarde een randpunt was, neem de waarde naast
    # het randpunt
    if index == 1
        index = 2
    elseif index == 7
        index = 6
    end

    # sla deze beginwaarde, afstand en periode op
    results = ["Xi" => [X[i] for i = (index-1):(index+1)],
               "Ei" => [d[i] for i = (index-1):(index+1)],
               "Ti" => [dists[i][2] for i = (index-1):(index+1)]];

    # test of de precisie nog niet om een kleinere stapgrootte vraagt
    if mini < delta_arr[h_index]
        h_index += 1;
        # we kunnen niet precieser binnen redelijke tijd, return de gevonden
        # waarden
        if h_index > num_h
            return results
        end
        # we kunnen precieser
        h = h_arr[h_index];
        H, n1, n2, N1, N2 = init_N(h, guess_T, 4);
    end

    # in de for-loop worden slechts 4 oplossingen opnieuw berekend,
    # pas de arrays voor pmap aan
    H = fill(h, 4);
    N1 = fill(n1, 4);
    N2 = fill(n2, 4);
    F = fill(f, 4);

    # verbeter de waarde voor y0 rounds keer
    for i=2:rounds
        # kopieer de waarden van de vorige iteratie naar hun nieuwe locatie
        X[1], X[4], X[7] = (X[index-1], X[index], X[index+1])
        dists[1], dists[4], dists[7] = (dists[index-1], dists[index],
                                        dists[index+1])
        # bereken de nieuwe X-waarden
        X[2], X[3] = ((2*X[1] + X[4])/3, (X[1] + 2*X[4])/3)
        X[5], X[6] = ((2*X[4] + X[7])/3, (X[4] + 2*X[7])/3)
        # bereken de nieuwe afstanden en periodes
        dists[2], dists[3], dists[5], dists[6] = pmap(RK4dists,
                F, [X[2], X[3], X[5], X[6]], H, N1, N2)

        # zoek de beste beginvoorwaarde en sla deze op, samen met de
        # afstand en periode
        d = [dist[1] for dist in dists];
        # zoek de index van de kleinste afstand
        mini = minimum(d);
        index = findfirst(isequal(mini), d);

        # als de beste waarde een randpunt was, neem de waarde naast
        # het randpunt
        if index == 1
            index = 2
        elseif index == 7
            index = 6
        end
        # sla deze beginwaarde, afstand en periode op
        push!(results,["Xi" => [X[i] for i = (index-1):(index+1)],
                       "Ei" => [d[i] for i = (index-1):(index+1)],
                       "Ti" => [dists[i][2] for i = (index-1):(index+1)]]);

        # test of de precisie nog niet om een kleinere stapgrootte vraagt
        if mini < delta_arr[h_index]
            h_index += 1;
            # we kunnen niet precieser binnen redelijke tijd, return de gevonden
            # waarden
            if h_index > num_h
                return results
            end
            # we kunnen precieser
            h = h_arr[h_index];
            H, n1, n2, N1, N2 = init_N(h, guess_T, 4);
        end
    end
    return results
end
