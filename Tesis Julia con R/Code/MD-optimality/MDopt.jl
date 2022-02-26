using CSV, DataFrames, NamedArrays, Combinatorics, LinearAlgebra, StatsBase
function MDopt(;X, y, Xcand, nMod, p_mod, fac_mod, nFDes, max_int, g, Iter, nStart, top)

    # To drop columns
    dropcol(M::NamedArray, j) = M[:, deleteat!(collect(axes(M, 2)), j)]
    dropcol(M::DataFrame, j) = M[:, deleteat!(collect(axes(M, 2)), j)]

    # Para quitar el primer elemento de un string
    snip(s::String) = s[nextind(s,1):end]
    snip(s::SubString) = s[nextind(s,1):end]

    # 1
    n = size(y, 1)
    fac = size(X, 2) - 1
    k = size(fac_mod, 2)

    Si = Vector{Float64}(undef, nMod)
    Xi = map(_ -> DataFrame(), 1:nMod) # Vector de Dataframes
    betai = Vector{Vector{Float64}}(undef,nMod)
    gammai = Matrix{Float64}[]  # Vector de matrices
    efectos = Vector{Vector{Int8}}(undef,nMod) # Vector de vectores

    models = zeros(nMod, fac)
    for i = 1:nMod
        aux = fac_mod[i, :]
        filter!(e -> e != 0, aux)
        models[i, aux] .= 1
    end

    Xfac = dropcol(X, 1)
    Xc = dropcol(Xcand, 1)

# 2
    models = NamedArray(models)
    if max_int > 1
        comb = hcat(collect(combinations(1:fac,2))...)
        mat = zeros(Int8, size(models, 1), size(comb, 2))
        for j = 1:size(comb, 2)
            fac1 = comb[1, j]
            fac2 = comb[2, j]
            aux = models[:, fac1] + models[:, fac2] .== 2
            filter!(e -> e != 0, aux)
            cols = names(aux, 1)
            cols = parse.(Int8, cols)
            mat[cols, j] .= 1
            Xfac = hcat(Xfac, Xfac[:, fac1].*Xfac[:, fac2], makeunique = true)
            Xc = hcat(Xc, Xc[:, fac1].*Xc[:, fac2])
        end
        mat = NamedArray(mat)
        models = hcat(models, mat)
    end

# 3
    if max_int > 2
        comb = hcat(collect(combinations(1:fac,3))...)
        mat = zeros(Int8, size(models, 1), size(comb, 2))
        for j = 1:size(comb, 2)
            fac1 = comb[1, j]
            fac2 = comb[2, j]
            fac3 = comb[3, j]
            aux = models[:, fac1] + models[:, fac2] + models[:, fac3] .== 3
            filter!(e -> e != 0, aux)
            cols = names(aux, 1)
            cols = parse.(Int8, cols)
            mat[cols, j] .= 1
            Xfac = hcat(Xfac, Xfac[:, fac1].*Xfac[:, fac2].*Xfac[:, fac3], makeunique = true)
            Xc = hcat(Xc, Xc[:, fac1].*Xc[:, fac2].*Xc[:, fac3])
        end
        mat = NamedArray(mat)
        models = hcat(models, mat)
    end

    Xfac = hcat(X[:, 1], Xfac, makeunique = true)
    Xc = hcat(Xcand[:, 1], Xc)
    models = hcat(ones(nMod), models)

    # 4
    models = NamedArray(models)
    for i=1:nMod
        aux = models[i, :] .== 1
        filter!(e -> e != 0, aux)
        cols = names(aux, 1)
        cols = parse.(Int8, cols)
        efectos[i] = cols
        tam = length(efectos[i])
        Xi[i] = insertcols!(Xfac[:, efectos[i]], 1,  :efectos => ones(n) , makeunique=true)

        mat = zeros(Int8, tam + 1, tam + 1)
        if size(mat, 1) > 1
            coord = hcat(2:size(mat, 1), 2:size(mat, 1))
            for k=1:(size(mat, 1) - 1)
                mat[coord[k, 1], coord[k, 2]] = 1
            end
        end

        push!(gammai, (1 / g^2)*mat)
        betai[i] = inv(gammai[i] + Matrix(Xi[i])' * Matrix(Xi[i])) * Matrix(Xi[i])' * y
        Si[i] = (y - Matrix(Xi[i]) * betai[i])' * (y - Matrix(Xi[i]) * betai[i]) + betai[i]'*gammai[i] * betai[i]
    end
   # 5
    function MDr(extra, X_mat)
        nex = length(extra)
        y_gorro_estrella = Vector{Vector{Float64}}(undef,nMod)
        V_estrella = Matrix{Float64}[]

        if typeof(X_mat) != Array
            X_mat = Array(X_mat)
        end

        for i=1:nMod
            Xiestrella = hcat(ones(nex), X_mat[extra, efectos[i]])
            y_gorro_estrella[i] = Xiestrella * betai[i]
            push!(V_estrella, diagm(0 => ones(nex)) + Xiestrella * inv(gammai[i] + Matrix(Xi[i])'*Matrix(Xi[i])) * Xiestrella')
        end

        MD = 0
        m = 1:nMod

        for i in m
            for j in m[m .!= i]
                MD = MD + p_mod[i]*p_mod[j]*(-nex +
                    sum(diag(inv(V_estrella[j])*V_estrella[i])) +
                    (n - 1)*(y_gorro_estrella[i] - y_gorro_estrella[j])'*
                    inv(V_estrella[j])*(y_gorro_estrella[i] - y_gorro_estrella[j])/Si[i])
            end
        end
        MD = MD * .5
        return(MD)
    end

    # 6
    df_MD = DataFrame(DesignPoints = [SubString("0")], MD = Float64(0))
    for j in 1:nStart
        extra = sample(1:size(Xcand, 1), nFDes, replace = true)
        iter = 1
        last_out = 0
        last_in = 1
        while last_out != last_in && iter < Iter
            dp = string(sort(extra))
            dp = snip(dp)  # Quita el primer parentesis
            dp = chop(dp)  # Quita el ultimo parentesis

            aux = df_MD[:, 1] in dp
            if aux != false
                break
            end

            df_MD = append!(df_MD, DataFrame(DesignPoints = [dp], MD = round(MDr(extra, Xc); digits = 2)))

            op = Vector{Float64}(undef, nFDes)
            for i in 1:nFDes
                op[i] = MDr(extra[1:end .!= i], Xc)
            end
            index = findmax(op)[2]
            last_out = extra[index]
            extra = extra[1:end .!= index]

            op = Vector{Float64}(undef, size(Xcand, 1))
            for i in 1:size(Xcand, 1)
                op[i] = MDr([extra; i],  Xc)
            end
            last_in = findmax(op)[2]
            append!(extra, last_in)

            iter = iter + 1
        end
    end

    df_MD = unique(df_MD, :DesignPoints)
    df_MD = sort(df_MD, [:MD], rev = true)
    df_MD = df_MD[1:top, :]

    # X
    aux_names = [Symbol("bkl")]
    append!(aux_names, [Symbol("F$i") for i in 1:fac])
    rename!(X, aux_names)

    # Xcand
    #aux_names = "F" .* string.(collect(1:1:fac))
    #pushfirst!(aux_names, "bkl")
    #setnames!(Xcand, aux_names, 2)

    # p_mod
    #aux_names = "M" .* string.(collect(1:1:nMod))
    #p_mod = NamedArray(p_mod)
    #setnames!(p_mod, aux_names, 1)

    # fac_mod
    #setnames!(fac_mod, aux_names, 1)

    # MDtop
    MDtop = df_MD[:, :MD]

    # DEStop
    DEStop = DataFrame()
    DEStop.DesignPoints = df_MD[:, :DesignPoints]

    return(df_MD)
end
