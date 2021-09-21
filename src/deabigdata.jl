# This file contains functions for the Radial Big Data Model, following Khezrimotlagh et al. (2019) algorithm (we use KZCT as acronym here)
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
    X_unselected = @view X[Not(dmus_selected), :]
    Y_unselected = @view Y[Not(dmus_selected), :]

    # create the pre-scores matrix (one column with the pre-score and the other to store the initial index)
    pre_scores = zeros(length(dmus_unselected))
    
    @inbounds for i in 1:size(X_unselected, 1)
        for x in 1:m 
            xquant = quantile!(X[:,x], 0.01:0.01:1)
            for k in 1:100 
                if X_unselected[i,x] <= xquant[k]
                    pre_scores[i] = pre_scores[i] + 1 
                end
            end
        end
        for y in 1:s 
            yquant = quantile!(Y[:,y], 0.01:0.01:1)
            for k in 1:100 
                if Y_unselected[i,y] >= yquant[k]
                    pre_scores[i] = pre_scores[i] + 1 
                end
            end
        end
    end

    # sort pre_scores to get the best DMUs in unselected DMUs 
    score_perm = sortperm(pre_scores, rev = true)
    dmus_unselected = dmus_unselected[score_perm]
    pre_scores = pre_scores[score_perm, :]

    # add the number of DMUs needed to fill out the subsample 
    z = p - length(dmus_selected)
    new_dmus = dmus_unselected[1:z]
    dmus_selected = vcat(dmus_selected, new_dmus)
    dmus_unselected = initial_index[Not(dmus_selected)]

    # Creates the subsets Dˢ and D_excluding_Dˢ
    return dmus_selected, dmus_unselected
end

"""
    bestpracticesfinder(Subset)

Identify best-practices in the subsample selected
"""
function bestpracticesfinder(scores::Vector, orient::Symbol, atol::Float64)
    index_bestpractices = Vector{Int64}()
    for i in 1:length(scores)
        if orient == :Input
            if scores[i] >= 1 - atol
                push!(index_bestpractices, i)
            end
        elseif orient == :Output
            if scores[i] <= 1 + atol
                push!(index_bestpractices, i)
            end
        end
    end
    return index_bestpractices
end

"""
    deabigdata(X, Y)

Compute the big data radial model using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `atol=1e-6`: tolerance for DMU to be considered efficient.
- `names`: a vector of strings with the names of the decision making units.
"""
function deabigdata(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
        orient::Symbol = :Input, rts::Symbol = :CRS, slack::Bool = true, atol::Float64 = 1e-6,
        names::Union{Vector{String},Nothing} = nothing,
        optimizer::Union{DEAOptimizer,Nothing} = nothing, progress::Bool = true)
    
    #=
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
    =#

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
    disposX = :Strong
    disposY = :Strong

    # Select initial subset Dˢ using Khezrimotlagh et al. (2019) algorithm 
    Dˢ, D_excluding_Dˢ = initialsubset(X, Y, n, s, m)

    # Find the best-practices Bˢ in Dˢ
    Dˢ_evaluation = dea(X[Dˢ, :], Y[Dˢ, :], orient = orient, rts = rts,slack = slack, 
                    disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)

    index_bestpractices = bestpracticesfinder(efficiency(Dˢ_evaluation), orient, atol)
    
    length(index_bestpractices) > 0 || throw(ErrorException("No efficient DMUs found in the initial subset. Consider increasing the tolerance `atol`"))

    Bˢ = Dˢ[index_bestpractices]

    # Find exterior DMUs in D excluding Dˢ respect to the hull of Bˢ 
    D_excluding_Dˢ_evaluation = dea(X[D_excluding_Dˢ, :], Y[D_excluding_Dˢ, :], orient = orient,
                        rts = rts, slack = slack, Xref = X[Bˢ, :], Yref = Y[Bˢ, :],
                        disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)

                    
    index_exteriors = bestpracticesfinder(efficiency(D_excluding_Dˢ_evaluation), orient, atol)     

    E = D_excluding_Dˢ[index_exteriors]

    # If E is not empty, then we have to find best practice DMUs in D^S Union E, otherwise, F = B_S and all DMUs are already evaluated 
    if size(E,1) > 0
        # Create the subset Dˢ_union_E        
        Dˢ_union_E = vcat(Dˢ, E)
        
        sort!(Dˢ_union_E, rev = false)

        Dˢ_union_E_evaluation = dea(X[Dˢ_union_E, :], Y[Dˢ_union_E, :], orient = orient, rts = rts,slack = slack,
                        disposX = disposX, disposY = disposY,optimizer = optimizer, progress = progress)
        
        # Find the index of best practices DMUs 
        index_bestpractices = bestpracticesfinder(efficiency(Dˢ_union_E_evaluation), orient, atol)
        
        # Best practices DMUs are the efficient DMUs F of the sample D
        F = Dˢ_union_E[index_bestpractices]

        # # We then evaluate the rest of DMUs 
        initial_index = collect(1:n)
        
        D_excluding_Dˢ_union_E = initial_index[Not(Dˢ_union_E)]

        D_excluding_Dˢ_union_E_evaluation = dea(X[D_excluding_Dˢ_union_E, :], Y[D_excluding_Dˢ_union_E, :], orient = orient,
                            rts = rts, slack = slack, Xref = X[F, :], Yref = Y[F, :],
                            disposX = disposX, disposY = disposY, optimizer = optimizer, progress = progress)

        subset1  = Dˢ_union_E
        results1 = Dˢ_union_E_evaluation
        subset2  = D_excluding_Dˢ_union_E
        results2 = D_excluding_Dˢ_union_E_evaluation   
    else
        F = Bˢ

        subset1  = Dˢ
        results1 = Dˢ_evaluation
        subset2  = D_excluding_Dˢ
        results2 = D_excluding_Dˢ_evaluation
    end

    # Get results
    resultsperm = sortperm(vcat(subset1, subset2), rev = false)
    effi = vcat(efficiency(results1), efficiency(results2))
    effi = effi[resultsperm]

    if slack
        slackX = vcat(slacks(results1, :X), slacks(results2, :X))
        slackY = vcat(slacks(results1, :Y), slacks(results2, :Y))

        slackX = slackX[resultsperm, :]
        slackY = slackY[resultsperm, :]
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    Xtarget = vcat(targets(results1, :X), targets(results2, :X))
    Ytarget = vcat(targets(results1, :Y), targets(results2, :Y))

    Xtarget = Xtarget[resultsperm, :]
    Ytarget = Ytarget[resultsperm, :]

    lambdas1 = peersmatrix(results1)
    lambdas2 = peersmatrix(results2)

    index_lambda_to_keep = findall(x -> x in F, subset1)
    lambdaeff = vcat(lambdas1[:,index_lambda_to_keep],lambdas2)
    lambdaindex = vcat(subset1, subset2)        
    lambda_matrix = lambdaeff[resultsperm, :]  

    # Create the sparse matrix for lambdas 
    lambdaeff = spzeros(n, n)
    lambdaeff[:, convert.(Int, lambdaindex[index_lambda_to_keep])] = lambda_matrix

    # return results
    return RadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)  
end


