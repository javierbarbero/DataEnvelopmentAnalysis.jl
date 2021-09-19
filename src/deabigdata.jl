# This file contains functions for the Khezrimotlagh et al. (2019) algorithm (we use KZCT as acronym here)
# D. Khezrimotlagh, J. Zhu and W.D. Cook et al. / European Journal of Operational Research 274 (2019) 1047–1054

"""
    initialsubset(X,Y)

    Creates the initial subset as in KZCT algorithm:
    1. Determines p the size of the subset and create an empty vector or matrix to fill out 
    2. Determines the minimum value for each input m and the maximum value for each output s 
    3. Select a subset with one DMU per input, with DMU's input equal to the minimum value for this input m and one DMU per output
    with DMU's output equal to the maximum value for this output s: m + s DMU's selected
    4. To fill out the subset up to p DMUs, implement the following algorithm:
        For each unselected DMU_j,
        {u = 0,
            {For each i and k_i,
                If x_ij <= k_ith percentile of x^i then
                    u = u + 1}
                {For each r and k'_r 
                    If y_rj >= k'th percentile of y^r then 
                        u = u + 1}
        pre_score_j = u }
    
        Where x^i and y^r represent the ith inputs and the rth outputs of all DMUS.
        The indexes k_i and k'_r are changed from 0 to 100 and refer to the percentiles of x^i and y^r respectively.
        Sorting DMUs in descending order by the assigned pre-scores, the remaining DMUs (to construct the subsample with size p)
        are selected as those having the greatest pre-scores. 

"""
function initialsubset(X::Union{Vector,Matrix},Y::Union{Vector,Matrix}, n::Int64, s::Int64, m::Int64)
    # Determines the size of the initial subset Dˢ
    p = convert(Int64, round(sqrt(n)))

    # Select DMU's with the minimum input or the maximum output
    dmus_selected = Vector{Int64}()
    for x in 1:m
        push!(dmus_selected, argmin(X[:,x]))
    end

    for y in 1:s
        push!(dmus_selected, argmax(Y[:,y]))
    end
    unique!(dmus_selected)

    initial_index = collect(1:n)
    dmus_unselected = initial_index[Not(dmus_selected)]

    # Algorithm implementation to fill out the subsample up to the size p 

    # select unselected DMUs 
    X_unselected = X[Not(dmus_selected), :]
    Y_unselected = Y[Not(dmus_selected), :]

    # create the pre-scores matrix (one column with the pre-score and the other to store the initial index)
    pre_scores = Vector(zeros(length(dmus_unselected)))
    pre_scores = hcat(pre_scores, dmus_unselected)
    
    for i in 1:size(X_unselected, 1)
        for x in 1:m 
            xquant = quantile!(X[:,x], 0.01:0.01:1)
            for k in 1:100 
                if X_unselected[i,x] <= xquant[k]
                    pre_scores[i, 1] = pre_scores[i, 1] + 1 
                end
            end
        end
        for y in 1:s 
            yquant = quantile!(Y[:,y], 0.01:0.01:1)
            for k in 1:100 
                if Y_unselected[i,y] >= yquant[k]
                    pre_scores[i, 1] = pre_scores[i, 1] + 1 
                end
            end
        end
    end

    # sort pre_scores to get the best DMUs in unselected DMUs 
    pre_scores = pre_scores[sortperm(pre_scores[:, 2], rev = true), :]

    # add the number of DMUs needed to fill out the subsample 
    z = p - length(dmus_selected)
    new_dmus = pre_scores[1:z,2]
    new_dmus = convert.(Int64, new_dmus)
    dmus_selected = vcat(dmus_selected, new_dmus)

    # Creates the subsets Dˢ and D_excluding_Dˢ
    Dˢ = [dmus_selected X[dmus_selected,:] Y[dmus_selected,:]]
    D_excluding_Dˢ = [initial_index[Not(dmus_selected)] X[Not(dmus_selected),:] Y[Not(dmus_selected),:]]
    return Dˢ, D_excluding_Dˢ
end

"""
    bestpracticesfinder(Subset)

    Identify best-practices in the subsample selected
"""
function bestpracticesfinder(scores::Vector, orient::Symbol)
    index_bestpractices = Vector{Int64}()
    for i in 1:length(scores)
        if orient == :Input
            if scores[i] >= 0.99
                push!(index_bestpractices, i)
            end
        elseif orient == :Output
            if scores[i] <= 1.01
                push!(index_bestpractices, i)
            end
        end
    end
    return index_bestpractices
end


"""
    getresults(Subset, evaluation)

    Add results to the matrix
"""
function getresults(subset::Matrix, evaluation::RadialDEAModel, n::Int64, m::Int64, s::Int64, slack::Bool)
    scores = efficiency(evaluation)
    if slack
        slacksX = slacks(evaluation, :X)
        slacksY = slacks(evaluation, :Y)
    else
        slacksX = zeros(evaluation.n, m)
        slacksY = zeros(evaluation.n, s)
    end
    lambdas = evaluation.lambda

    subset = [subset[:,1:m+1+s] scores slacksX slacksY]
    return subset, lambdas
end


"""
    deabigdata(X, Y)

    Apply Radial DEA Model using KZCT algorithm
    The algorithm is such as:
    1. Start 
    2. D <- get a sample of DMUs 
    3. D^S <- select a subsample of D 
    4. B^S <- find best-practices in D^S 
    5. E <- find exterior DMUs in D excluding D^S respect to the hull of B^S 
    6. If E = {}, then F = B^S and all DMUs are altready evaluated (go to step 7)
    Otherwise 
        6.1. F <- find best-practice DMUs in D^S Union E 
        6.2. Evaluate DMUs in D excluding (D^S Union E) respoect to the F's hull 
    7. End 

    Where F is the set of efficient DMUs. 

    # Optional Arguments
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `disposX=:Strong`: chooses strong disposability of inputs. For weak disposability choose `:Weak`.
- `disposY=:Strong`: chooses strong disposability of outputs. For weak disposability choose `:Weak`.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deabigdata(X, Y)
Radial DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
────────────────────────────────────────────────────────
    efficiency       slackX1       slackX2       slackY1
────────────────────────────────────────────────────────
1     1.0       -1.10378e-8   -1.55523e-8   -2.14142e-8
2     0.62229   -3.8778e-11   -3.01925e-11  -3.79858e-11
3     0.819856   4.18927e-11   3.39628e-11   4.72502e-11
4     1.0        3.22358e-14   8.91749e-14   1.31448e-13
5     0.310371   2.29024e-11   3.03865e-11   3.26362e-11
6     0.555556   4.44444      -6.72088e-12  -2.42811e-12
7     1.0        6.40683e-11  -9.57711e-13   5.13434e-12
8     0.757669   2.82994e-10  -2.42352e-17  -6.78179e-17
9     0.820106   1.64021      -1.38876e-7   -1.86896e-7
10    0.490566   2.66381e-11   7.62939e-12   1.75776e-11
11    1.0       -1.47048e-11   4.0          -2.66588e-11
────────────────────────────────────────────────────────
```
"""

function deabigdata(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
        orient::Symbol = :Input, rts::Symbol = :CRS, slack::Bool = true,
        disposX::Symbol = :Strong, disposY::Symbol = :Strong,
        names::Union{Vector{String},Nothing} = nothing,
        optimizer::Union{DEAOptimizer,Nothing} = nothing, progress::Bool = true)
    
    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    if nx != ny
        throw(DimensionMismatch("number of@rows in X and Y ($nx, $ny) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    n = nx

#################################### SELECT A SUBSET Dˢ ########################################################    
    Dˢ, D_excluding_Dˢ = initialsubset(X, Y, n, s, m)

    # Fix the coordinates in the matrix 
    coordX = [2, 2+m-1]
    coordY = [2+m, 2+m+s-1]
    coordscores = [2+m+s]
    coordslackX = [coordscores[1]+1,coordscores[1]+1+m-1]
    coordslackY = [coordslackX[2]+1, coordslackX[2]+1+s-1]
#################################### FIND BEST-PRACTICES Bˢ ########################################################    

    # Find the best-practices Bˢ in Dˢ
    evaluation = dea(Dˢ[:,coordX[1]:coordX[2]], Dˢ[:,coordY[1]:coordY[2]], orient = orient, rts = rts,slack = slack, 
                    disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)
    Dˢ, lambdas_Dˢ = getresults(Dˢ, evaluation, n, m, s, slack)

    index_bestpractices = bestpracticesfinder(Dˢ[:,coordscores[1]], orient)

    Bˢ = Dˢ[index_bestpractices,:]

#################################### FIND EXTERIORS E ########################################################    

    # Find exterior DMUs in D excluding Dˢ respect to the hull of Bˢ 
    evaluation = dea(D_excluding_Dˢ[:,coordX[1]:coordX[2]], D_excluding_Dˢ[:,coordY[1]:coordY[2]], orient = orient,
                        rts = rts, slack = slack, Xref = Bˢ[:,coordX[1]:coordX[2]], Yref = Bˢ[:,coordY[1]:coordY[2]],
                        disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)

    D_excluding_Dˢ, lambdas_D_excluding_Dˢ = getresults(D_excluding_Dˢ, evaluation, n, m, s, slack)
                    
    index_exteriors = bestpracticesfinder(D_excluding_Dˢ[:,coordscores[1]], orient)     

    E = D_excluding_Dˢ[index_exteriors, :]

#################################### IF EXTERIORS E, FIND D IN Dˢ_union_E ########################################################    

    # If E is not empty, then we have to find best practice DMUs in D^S Union E, otherwise, F = B_S and all DMUs are already evaluated 
    if size(E,1)>0
        # Create the subset Dˢ_union_E
        Dˢ_union_E = vcat(Dˢ, E)
        Dˢ_union_E = Dˢ_union_E[sortperm(Dˢ_union_E[:, 1], rev = false), :]

        evaluation = dea(Dˢ_union_E[:,coordX[1]:coordX[2]], Dˢ_union_E[:,coordY[1]:coordY[2]], orient = orient, rts = rts,slack = slack,
                        disposX = disposX, disposY = disposY,optimizer = optimizer, progress = progress)
        
        Dˢ_union_E, lambdas_Dˢ_union_E = getresults(Dˢ_union_E, evaluation, n, m, s, slack)

        # Find the index of best practices DMUs 
        index_bestpractices = bestpracticesfinder(Dˢ_union_E[:,coordscores[1]], orient)
        
        # Best practices DMUs are the efficient DMUs F of the sample D
        F = Dˢ_union_E[index_bestpractices,:]

        # # We then evaluate the rest of DMUs 
        initial_index = collect(1:n)
        
        D_excluding_Dˢ_union_E = [initial_index[Not(Dˢ_union_E[:,1])] X[Not(Dˢ_union_E[:,1]),:] Y[Not(Dˢ_union_E[:,1]),:]]
        
        evaluation = dea(D_excluding_Dˢ_union_E[:,coordX[1]:coordX[2]], D_excluding_Dˢ_union_E[:,coordY[1]:coordY[2]], orient = orient,
                            rts = rts, slack = slack, Xref = F[:,coordX[1]:coordX[2]], Yref = F[:,coordY[1]:coordY[2]],
                            disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)

        D_excluding_Dˢ_union_E, lambdas_D_excluding_Dˢ_union_E = getresults(D_excluding_Dˢ_union_E, evaluation, n, m, s, slack)

        results = vcat(Dˢ_union_E, D_excluding_Dˢ_union_E)
        results = results[sortperm(results[:, 1], rev = false), :]

        index_lambda_to_keep = findall(x -> x in F[:,1],Dˢ_union_E[:,1])
        lambdaeff = vcat(lambdas_Dˢ_union_E[:,index_lambda_to_keep],lambdas_D_excluding_Dˢ_union_E)
        lambdaindex = vcat(Dˢ_union_E[:,1], D_excluding_Dˢ_union_E[:,1])
        lambda_matrix = lambdaeff[sortperm(lambdaindex, rev = false), :]
    else
        F = Bˢ
        results = vcat(Dˢ, D_excluding_Dˢ)
        results = results[sortperm(results[:, 1], rev = false), :]

        index_lambda_to_keep = findall(x -> x in F[:,1],D_excluding_Dˢ[:,1])
        lambdaeff = vcat(lambdas_Dˢ[:,index_lambda_to_keep],lambdas_D_excluding_Dˢ)
        lambdaindex = vcat(lambdas_Dˢ[:,1], lambdas_D_excluding_Dˢ[:,1])
        lambda_matrix = lambdaeff[sortperm(lambdaindex, rev = false), :]
    end

    # Create the sparse matrix for lambdas 
    lambdaeff = spzeros(n, n)
    lambdaeff[:, convert.(Int, lambdaindex[index_lambda_to_keep])] = lambda_matrix

    effi = results[:,coordscores[1]]
    if slack
        slackX = results[:, coordslackX[1]:coordslackX[2]]
        slackY = results[:, coordslackY[1]:coordslackY[2]]
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end    

    # Get X and Y targets
    if orient == :Input
        Xtarget = X .* effi
        Ytarget = Y
    elseif orient == :Output
        Xtarget = X
        Ytarget = Y .* effi
    end

    if slack        
        Xtarget = Xtarget - slackX
        Ytarget = Ytarget + slackY
    else
        if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
        if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end
    end

    # return results
    return RadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)  
end


