# This file contains functions for the Khezimotlagh et al. (2019) algorithm (we use KZCT as acronym here)
# D. Khezrimotlagh, J. Zhu and W.D. Cook et al. / European Journal of Operational Research 274 (2019) 1047–1054
"""
    Subset
Model subset as implemented in KZCT algorithm (2019)
"""
Base.@kwdef mutable struct Subset 
    X::Union{Matrix,Vector} # Inputs subset of selected DMUs
    Y::Union{Matrix,Vector} # Outputs subset of selected DMUs 
    indexDMU::Vector{Int64} # Initial index of selected DMUs
end

"""
    initialsubset(X, Y)

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
function initialsubset(X::Union{Vector,Matrix}, Y::Union{Vector,Matrix})
    
    # Step 1 
    # Number of DMUs
    n = size(X,1)
    
    # Size of the subset to create 
    p = convert(Int64, round(sqrt(n)))

    # number of inputs 
    m = size(X,2)

    # number of ouputs 
    s = size(Y,2)
    
    # Step 2 
    # Determine the minimum value for each input m 
    X_min = Vector(zeros(m))
    for i in 1:m 
        X_min[i] = minimum(X[:,i])
    end
    # Determine the maximum value for each output s 
    Y_max = Vector(zeros(m))
    for i in 1:s 
        Y_max[i] = maximum(Y[:,i])
    end

    # Step 3 
    # Select a subset with m + s DMUs 
    dmus_selected = Vector{Int64}()
    for x in 1:m 
        counter_x = 0 
        for i in 1:n 
            if (X[i,x] <= X_min[x] && counter_x == 0)
                push!(dmus_selected, i)
                counter_x = counter_x + 1 
            else 
                nothing 
            end
        end
    end
    for y in 1:s
        counter_y = 0 
        for i in 1:n 
            if (Y[i, y] >= Y_max[y] && counter_y == 0)
                push!(dmus_selected, i)
                counter_y = counter_y + 1 
            else
                nothing 
            end
        end
    end

    unique!(dmus_selected)
    initial_index = [i for i in 1:n]
    dmus_unselected = initial_index[Not(dmus_selected)]

    # Step 4
    # Algorithm implementation to fill out the subsample up to the size p 

    # select unselected DMUs 
    X_unselected = X[Not(dmus_selected), :]
    Y_unselected = Y[Not(dmus_selected), :]

    # create the pre-scores matrix (one column with the pre-score and the other to store the initial index)
    pre_scores = Vector(zeros(length(dmus_unselected)))
    pre_scores = hcat(pre_scores, dmus_unselected)
    
    for i in 1:size(X_unselected, 1)
        for x in 1:m 
            for k in 1:100 
                if X_unselected[i,x] <= quantile!(X[:,x], k / 100)
                    pre_scores[i, 1] = pre_scores[i, 1] + 1 
                else
                    nothing
                end
            end
        end
        for y in 1:s 
            for k in 1:100 
                if Y_unselected[i,y] >= quantile!(Y[:,y], k / 100)
                    pre_scores[i, 1] = pre_scores[i, 1] + 1 
                else
                    nothing
                end
            end
        end
    end

    # sort pre_scores to get the best DMUs in unselected DMUs 
    pre_scores = pre_scores[sortperm(pre_scores[:, 2], rev = true), :]

    # add the number of DMUs needed to fill out the subsample 
    z = p - length(dmus_selected)
    new_dmus = pre_scores[1:z,2]
    dmus_selected = vcat(dmus_selected, new_dmus)

    D_S = Subset(X = X[dmus_selected,:], Y = Y[dmus_selected, :], indexDMU = dmus_selected)
    Non_D_S = Subset(X = X_unselected, Y = Y_unselected, indexDMU = dmus_unselected)

    return D_S, Non_D_S
end


"""
    bestpracticesfinder(Subset_1, orient, rts, slack)

    Identify best-practices in the subsample selected
"""

function bestpracticesfinder(D_S::Subset, orient::Symbol, rts::Symbol, slack::Bool)
    evaluation = deamodel(D_S.X, D_S.Y, orient = orient, rts = rts, slack = slack)

    scores = efficiency(evaluation)

    index_bestpractices = Vector{Int64}()
    for i in 1:length(scores)
        if scores[i] >= 1
            push!(index_bestpractices, i)
        else
            nothing
        end
    end

    best_practices = Subset(X = D_S.X[index_bestpractices,:], 
                            Y = D_S.Y[index_bestpractices,:], 
                            indexDMU = D_S.indexDMU[index_bestpractices,:])

    return best_practices, evaluation 
end


"""
    exteriorfinder(Subset_1,Subset2, orient, rts, slack)

    Find exterior DMUs in D \ D^S respect to the hull of B^S
"""

function exteriorsfinder(B_S::Subset, Non_D_S::Subset, orient::Symbol, rts::Symbol, slack::Bool)
    eval = deamodel(Non_D_S.X, Non_D_S.Y, Xref = B_S.X, Y_ref = B_S.Y, orient = orient, rts = rts, slack = slack)
    scores = efficiency(eval)

    index_exteriors = Vector{Int64}()
    for i in 1:length(scores)
        if scores[i] >= 1
            push!(index_exteriors, i)
        else
            nothing
        end
    end

    E = Subset(X = Non_D_S.X[index_exteriors,:], 
                            Y = Non_D_S.Y[index_exteriors,:], 
                            indexDMU = Non_D_S.indexDMU[index_exteriors,:])
    
    return E, eval 

end

"""
    reformatlambda(results_1, results_2, index_1, index_2)

    Remormat lambdas from results to return lambdas in the same order / format than the initial sample D 
"""
function reformatlambda(results_Non_D_S_Union_E::RadialDEAModel, results_D_S_Union_E::RadialDEAModel, F::Subset, D_S_Union_E::Subset)

    lambda_1 = results_Non_D_S_Union_E.lambda # F is the hull 
    lambda_2 = results_D_S_Union_E.lambda # F + other are the hull here, we have to drop those which are not hull   
    
    index_to_keep = findall(x -> x in F.indexDMU, results_D_S_Union_E.indexDMU)
    lambda_1 = lambda_1[:,index_to_keep]

    lambda = vcat(lambda_1,lambda_2)

    return lambda

end

"""
    bigdata(X, Y)

    Apply Radial DEA Model using KZCT algorithm
    The algorithm is such as:
    1. Start 
    2. D <- get a sample of DMUs 
    3. D^S <- select a subsample of D 
    4. B^S <- find best-practices in D^S 
    5. E <- find exterior DMUs in D \ D^S respect to the hull of B^S 
    6. If E = {}, then F = B^S and all DMUs are altready evaluated (go to step 7)
    Otherwise 
        6.1. F <- find best-practice DMUs in D^S Union E 
        6.2. Evaluate DMUs in D \ (D^S Union E) respoect to the F's hull 
    7. End 

    Where F is the set of efficient DMUs. 

    # Optional Arguments
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `disposX=:Strong`: chooses strong disposability of inputs. For weak disposability choose `:Weak`.
- `disposY=:Strong`: chooses strong disposability of outputs. For weak disposability choose `:Weak`.
"""

function bigdata(X::Union{Vector, Matrix}, Y::Union{Vector, Matrix}, orient::Symbol, rts::Symbol, slack::Bool = false)
    # Create and fill the subsample
    D_S, Non_D_S = initialsubset(X, Y)
    entire_index_DMU = vcat(D_S.indexDMU, Non_D_S.indexDMU)
    # Find the best-practices B_S in D_S
    B_S, results_D_S = bestpracticesfinder(D_S, orient, rts, slack)
    # Find exterior DMUs in D \ D^S respect to the hull of B^S 
    E, results_non_D_S = exteriorsfinder(B_S, Non_D_S, orient, rts, slack)
    # If E is not empty, then we have to find best practice DMUs in D^S Union E, otherwise, F = B_S and all DMUs are already evaluated 
    if size(E.X,1)>1
        D_S_Union_E = Subset(X = vcat(D_S.X, E.X), Y = vcat(D_S.X, E.X), indexDMU = vcat(D_S.indexDMU, E.indexDMU))
        F, results_D_S_Union_E = bestpracticesfinder(D_S_Union_E, orient, rts, slack)
        Non_D_S_Union_E = Subset(X = X[Not(D_S_Union_E.indexDMU),:], Y = Y[Not(D_S_Union_E.indexDMU),:], indexDMU = entire_index_DMU[Not(D_S_Union_E.indexDMU)])
        Other, results_Non_D_S_Union_E = exteriorsfinder(F, Non_D_S_Union_E, orient, rts, slack) 
        Non_F = Subset(X = X[Not(F.indexDMU),:], Y = Y[Not(F.indexDMU), :], indexDMU = entire_index_DMU[Not(F.indexDMU)])
        # reformat to get in the same order for effi
        effi_scores = vcat(efficiency(results_D_S_Union_E),efficiency(results_Non_D_S_Union_E))
        effi_index = vcat(D_S_Union_E.indexDMU, Non_D_S_Union_E.indexDMU)
        effi_matrix = hcat(effi_scores, effi_index)
        effi_matrix = effi_matrix[sortperm(effi_matrix[:, 2], rev = true), :]
        effi = effi_matrix[:,1]

        # reformat lambda
        lambdaeff = reformatlambda(results_Non_D_S_Union_E, results_D_S_Union_E, F, D_S_Union_E)

    else
        F = B_S
        Non_F = Subset(X = X[Not(F.indexDMU),:], Y = Y[Not(F.indexDMU), :], indexDMU = entire_index_DMU[Not(F.indexDMU)])
        # reformat to get in the same order for effi
        effi_scores = vcat(efficiency(results_D_S),efficiency(results_non_D_S))
        effi_index = vcat(D_S.indexDMU, Non_D_S.indexDMU)
        effi_matrix = hcat(effi_scores, effi_index)
        effi_matrix = effi_matrix[sortperm(effi_matrix[:, 2], rev = true), :]
        effi = effi_matrix[:,1]

        # reformat lambda
        lambdaeff = reformatlambda(results_non_D_S, results_D_S, B_S, D_S)
    end 
    
    return effi, lambdaeff 
end
