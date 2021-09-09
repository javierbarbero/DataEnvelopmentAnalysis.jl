# This file contains functions for the Khezrimotlagh et al. (2019) algorithm (we use KZCT as acronym here)
# D. Khezrimotlagh, J. Zhu and W.D. Cook et al. / European Journal of Operational Research 274 (2019) 1047–1054
"""
    Subset
A data structure representing a subset of DMUs.
"""
Base.@kwdef mutable struct Subset 
    X::Union{Matrix,Vector} = [] # Inputs subset of selected DMUs
    Y::Union{Matrix,Vector} = [] # Outputs subset of selected DMUs 
    indexDMU::Vector{Int64} = [] # Initial index of selected DMUs
    eff::Vector{Float64} = [] # Efficiency scores
    lambdaeff::SparseMatrixCSC{Float64, Int64} = spzeros(2,2) # Optimal weights
    slackX::Matrix = Matrix(zeros(2,2))
    slackY::Matrix = Matrix(zeros(2,2))
    Xtarget::Matrix = Matrix(zeros(2,2))
    Ytarget::Matrix = Matrix(zeros(2,2))
end

"""
    KZCTAlgorithm
A data structure representing the Khezrimotlagh et al. (2019) algorithm.
"""
Base.@kwdef mutable struct KZCTAlgorithm
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    disposX::Symbol
    disposY::Symbol
    slack::Bool
    optimizer::DEAOptimizer
    dmunames::Union{Vector{String},Nothing}
    eff::Vector = Vector{Float64}()
    slackX::Matrix = Matrix(zeros(2,2))
    slackY::Matrix = Matrix(zeros(2,2))
    lambda::SparseMatrixCSC{Float64, Int64} = spzeros(2,2)
    Xtarget::Matrix = Matrix(zeros(2,2))
    Ytarget::Matrix= Matrix(zeros(2,2))
    D::Subset = Subset()# Initial Sample  
    Dˢ::Subset = Subset() # Selected Sample 
    D_exluding_Dˢ::Subset = Subset() # Unselected sample 
    Bˢ::Subset = Subset() # Best-practices among Dˢ
    Ε::Subset = Subset() # Find exterior in D_exluding_Dˢ respect to the hull of Bˢ
    Dˢ_union_E::Subset = Subset() # Union Dˢ subset and E subset
    F::Subset = Subset() # Efficient DMUs 
    D_exluding_Dˢ_union_Ε::Subset = Subset() # Inefficient DMUs to be evaluated respect to the F's hull 
end


"""
    initialsubset!(KZCT_algorithm)

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
function initialsubset!(KZCT_algorithm::KZCTAlgorithm)
    
    X = KZCT_algorithm.D.X 
    Y = KZCT_algorithm.D.Y

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
    initial_index = KZCT_algorithm.D.indexDMU
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
    new_dmus = convert.(Int64, new_dmus)
    dmus_selected = vcat(dmus_selected, new_dmus)

    KZCT_algorithm.Dˢ = Subset(X = X[dmus_selected,:], Y = Y[dmus_selected, :], indexDMU = dmus_selected)
    KZCT_algorithm.D_exluding_Dˢ = excludingsubset(X, Y, dmus_selected)
    
end


"""
    bestpracticesfinder!(KZCT_algorithm, subset)

    Identify best-practices in the subsample selected
"""

function bestpracticesfinder!(KZCT_algorithm::KZCTAlgorithm, subset::Symbol)

    if subset == :Dˢ
        evaluation = dea(KZCT_algorithm.Dˢ.X, KZCT_algorithm.Dˢ.Y, 
                        orient = KZCT_algorithm.orient, rts = KZCT_algorithm.rts, 
                        slack = KZCT_algorithm.slack, optimizer = KZCT_algorithm.optimizer)
    elseif subset == :Dˢ_union_E
        evaluation = dea(KZCT_algorithm.Dˢ_union_E.X, KZCT_algorithm.Dˢ_union_E.Y, 
                        orient = KZCT_algorithm.orient, rts = KZCT_algorithm.rts, 
                        slack = KZCT_algorithm.slack, optimizer = KZCT_algorithm.optimizer)
    else 
        throw(ArgumentError("`subset` must be :Dˢ or :Dˢ_union_E"));
    end

    scores = efficiency(evaluation)

    index_bestpractices = Vector{Int64}()
    for i in 1:length(scores)
        if scores[i] >= 0.99
            push!(index_bestpractices, i)
        else
            nothing
        end
    end

    if subset == :Dˢ
        index_to_keep = KZCT_algorithm.Dˢ.indexDMU[index_bestpractices]
        KZCT_algorithm.Bˢ = createsubset(KZCT_algorithm.D.X, KZCT_algorithm.D.Y, index_to_keep)
        KZCT_algorithm.Bˢ.eff = scores[index_bestpractices]
        KZCT_algorithm.Dˢ.eff = scores 
        KZCT_algorithm.Bˢ.lambdaeff = evaluation.lambda 
        KZCT_algorithm.Dˢ.lambdaeff = evaluation.lambda
        KZCT_algorithm.Dˢ.slackX = slacks(evaluation, :X)
        KZCT_algorithm.Dˢ.slackY = slacks(evaluation, :Y)
        KZCT_algorithm.Dˢ.Xtarget = targets(evaluation, :X)
        KZCT_algorithm.Dˢ.Ytarget = targets(evaluation, :Y)
    else 
        index_to_keep = KZCT_algorithm.Dˢ_union_E.indexDMU[index_bestpractices]
        KZCT_algorithm.F = createsubset(KZCT_algorithm.D.X, KZCT_algorithm.D.Y, index_to_keep)
        KZCT_algorithm.F.eff = scores[index_bestpractices]
        KZCT_algorithm.Dˢ_union_E.eff = scores 
        KZCT_algorithm.F.lambdaeff = evaluation.lambda 
        KZCT_algorithm.Dˢ_union_E.lambdaeff = evaluation.lambda
        KZCT_algorithm.Dˢ_union_E.slackX = slacks(evaluation, :X)
        KZCT_algorithm.Dˢ_union_E.slackY = slacks(evaluation, :Y)
        KZCT_algorithm.Dˢ_union_E.Xtarget = targets(evaluation, :X)
        KZCT_algorithm.Dˢ_union_E.Ytarget = targets(evaluation, :Y)
    end
    
end


"""
    exteriorsfinder!(KZCT_algorithm, subset)

    Find exterior DMUs respect to the hull
"""

function exteriorsfinder!(KZCT_algorithm::KZCTAlgorithm, subset::Symbol)
    X = KZCT_algorithm.D.X 
    Y = KZCT_algorithm.D.Y 

    if subset == :Bˢ
        eval = dea(KZCT_algorithm.D_exluding_Dˢ.X, KZCT_algorithm.D_exluding_Dˢ.Y, 
                    orient = KZCT_algorithm.orient, rts = KZCT_algorithm.rts, slack = KZCT_algorithm.slack, 
                    Xref = KZCT_algorithm.Bˢ.X, Yref = KZCT_algorithm.Bˢ.Y, optimizer = KZCT_algorithm.optimizer)
    elseif subset == :F 
        eval = dea(KZCT_algorithm.D_exluding_Dˢ_union_Ε.X, KZCT_algorithm.D_exluding_Dˢ_union_Ε.Y, 
        orient = KZCT_algorithm.orient, rts = KZCT_algorithm.rts, slack = KZCT_algorithm.slack, 
        Xref = KZCT_algorithm.F.X, Yref = KZCT_algorithm.F.Y, optimizer = KZCT_algorithm.optimizer)
    else 
        throw(ArgumentError("`subset` must be :Bˢ or :F"));
    end

    scores = efficiency(eval)

    index_exteriors = Vector{Int64}()
    for i in 1:length(scores)
        if scores[i] >= 0.99
            push!(index_exteriors, i)
        else
            nothing
        end
    end

    if subset == :Bˢ
        index_to_keep = KZCT_algorithm.D_exluding_Dˢ.indexDMU[index_exteriors]
        KZCT_algorithm.Ε = createsubset(X, Y, index_to_keep)
        KZCT_algorithm.Ε.eff = scores[index_exteriors]
        KZCT_algorithm.Ε.lambdaeff = eval.lambda[index_exteriors]
        KZCT_algorithm.D_exluding_Dˢ.eff = scores 
        KZCT_algorithm.D_exluding_Dˢ.lambdaeff = eval.lambda 
        KZCT_algorithm.D_exluding_Dˢ.slackX = slacks(eval, :X)
        KZCT_algorithm.D_exluding_Dˢ.slackY = slacks(eval, :Y)
        KZCT_algorithm.D_exluding_Dˢ.Xtarget = targets(eval, :X)
        KZCT_algorithm.D_exluding_Dˢ.Ytarget = targets(eval, :Y)
    else
        index_to_keep = KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU[index_exteriors]
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.eff = scores 
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.lambdaeff = eval.lambda 
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.slackX = slacks(eval, :X)
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.slackY = slacks(eval, :Y)
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.Xtarget = targets(eval, :X)
        KZCT_algorithm.D_exluding_Dˢ_union_Ε.Ytarget = targets(eval, :Y)
    end

    

end



"""
    excludingsubset(X, Y, index_to_exclude)

    Creates the subset excluding the index 
"""

function excludingsubset(X::Union{Vector,Matrix}, Y::Union{Vector,Matrix}, index_to_exclude::Vector{Int64})
    entire_index = [i for i in 1:size(X,1)]

    index_to_keep = findall(x -> !(x in index_to_exclude), entire_index)

    return Subset(X = X[index_to_keep, :], Y = Y[index_to_keep,:], indexDMU = entire_index[index_to_keep])
end

"""
    unionsubset(X, Y, indexsubsetone, indexsubsettwo)

    Creates the subset unioning subset one and subset two
"""

function unionsubset(X::Union{Vector,Matrix}, Y::Union{Vector,Matrix}, indexsubsetone::Vector{Int64}, indexsubsettwo)
    union_index = vcat(indexsubsetone,indexsubsettwo)

    entire_index = [i for i in 1:size(X,1)]

    return Subset(X = X[union_index, :], Y = Y[union_index,:], indexDMU = entire_index[union_index])

end

"""
    createsubset(X,Y, index_to_include)

    Creates the subset based on original X and Y sample and the index of DMUs to keep
"""
function createsubset(X::Union{Vector,Matrix}, Y::Union{Vector,Matrix}, indexsubset::Vector{Int64})
    sort!(indexsubset)

    return Subset(X = X[indexsubset,:], Y = Y[indexsubset,:], indexDMU = indexsubset)
end
"""
    getscores!(KZCT_algorithm)

    Get the scores in the proper order of the sample D 
"""

function getscores!(KZCT_algorithm::KZCTAlgorithm,  Exterior_presence::Bool)
    if Exterior_presence == true 
        eff = vcat(KZCT_algorithm.Dˢ_union_E.eff, KZCT_algorithm.D_exluding_Dˢ_union_Ε.eff)
        eff_index = vcat(KZCT_algorithm.Dˢ_union_E.indexDMU, KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU)
        effi_matrix = hcat(eff, eff_index)
        effi_matrix = effi_matrix[sortperm(effi_matrix[:, 2], rev = false), :]
        KZCT_algorithm.eff = effi_matrix[:,1]
    else 
        eff = vcat(KZCT_algorithm.Dˢ.eff, KZCT_algorithm.D_exluding_Dˢ.eff)
        eff_index = vcat(KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
        effi_matrix = hcat(eff, eff_index)
        effi_matrix = effi_matrix[sortperm(effi_matrix[:, 2], rev = false), :]
        KZCT_algorithm.eff = effi_matrix[:,1]
    end

end

"""
    getlambdas!(KZCT_algorithm, Exterior_presence)

    Get the lambdas in the proper order of the sample D 
"""

function getlambdas!(KZCT_algorithm::KZCTAlgorithm,  Exterior_presence::Bool)
    if Exterior_presence == true 
        lambda_1 = KZCT_algorithm.Dˢ_union_E.lambdaeff # F is the hull 
        lambda_2 = KZCT_algorithm.D_exluding_Dˢ_union_Ε.lambdaeff # F + other are the hull here, we have to drop those which are not hull   
        index_to_keep = findall(x -> x in KZCT_algorithm.F.indexDMU, KZCT_algorithm.Dˢ_union_E.indexDMU)
        lambda_1 = lambda_1[:,index_to_keep]
        KZCT_algorithm.lambda = vcat(lambda_1,lambda_2)
    else 
        lambda_1 = KZCT_algorithm.Dˢ.lambdaeff # F is the hull 
        lambda_2 = KZCT_algorithm.D_exluding_Dˢ.lambdaeff # F + other are the hull here, we have to drop those which are not hull   
        index_to_keep = findall(x -> x in KZCT_algorithm.F.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
        lambda_1 = lambda_1[:,index_to_keep]
        KZCT_algorithm.lambda = vcat(lambda_1,lambda_2)
    end
end

"""
    getslacks!(KZCT_algorithm, slacks, Exterior_presence)

    Get the slacks in the proper order of the sample D 
"""

function getslacks!(KZCT_algorithm::KZCTAlgorithm, slacks::Symbol, Exterior_presence::Bool)

    if Exterior_presence == true 
        if slacks == :X
            slackX = vcat(KZCT_algorithm.Dˢ_union_E.slackX, KZCT_algorithm.D_exluding_Dˢ_union_Ε.slackX)
            slackX_index = vcat(KZCT_algorithm.Dˢ_union_E.indexDMU, KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU)
            slackX_matrix = hcat(slackX, slackX_index)
            slackX_matrix = slackX_matrix[sortperm(slackX_matrix[:, size(slackX,2)+1], rev = false), :]
            KZCT_algorithm.slackX = slackX_matrix[:,1:size(slackX,2)]

        elseif slacks == :Y
            slackY = vcat(KZCT_algorithm.Dˢ_union_E.slackY, KZCT_algorithm.D_exluding_Dˢ_union_Ε.slackY)
            slackY_index = vcat(KZCT_algorithm.Dˢ_union_E.indexDMU, KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU)
            slackY_matrix = hcat(slackY, slackY_index)
            slackY_matrix = slackY_matrix[sortperm(slackY_matrix[:, size(slackY,2)+1], rev = false), :]
            KZCT_algorithm.slackY = slackY_matrix[:,1:size(slackY,2)]
        else 
            throw(ArgumentError("`slacks` must be :X or :Y"));
        end
    else 
        if slacks == :X
            slackX = vcat(KZCT_algorithm.Dˢ.slackX, KZCT_algorithm.D_exluding_Dˢ.slackX)
            slackX_index = vcat(KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
            slackX_matrix = hcat(slackX, slackX_index)
            slackX_matrix = slackX_matrix[sortperm(slackX_matrix[:, size(slackX,2)+1], rev = false), :]
            KZCT_algorithm.slackX = slackX_matrix[:,1:size(slackX,2)]

        elseif slacks == :Y
            slackY = vcat(KZCT_algorithm.Dˢ.slackY, KZCT_algorithm.D_exluding_Dˢ.slackY)
            slackY_index = vcat(KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
            slackY_matrix = hcat(slackY, slackY_index)
            slackY_matrix = slackY_matrix[sortperm(slackY_matrix[:, size(slackY,2)+1], rev = false), :]
            KZCT_algorithm.slackY = slackY_matrix[:,1:size(slackY,2)]
        else 
            throw(ArgumentError("`slacks` must be :X or :Y"));
        end
    end
end

"""
    gettargets!(KZCT_algorithm, targets, Exterior_presence)

    Get the targets in the proper order of the sample D 
"""

function gettargets!(KZCT_algorithm::KZCTAlgorithm, targets::Symbol, Exterior_presence::Bool)
    if Exterior_presence == true 
        if targets == :X
            Xtarget = vcat(KZCT_algorithm.Dˢ_union_E.Xtarget, KZCT_algorithm.D_exluding_Dˢ_union_Ε.Xtarget)
            Xtarget_index = vcat(KZCT_algorithm.Dˢ_union_E.indexDMU, KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU)
            Xtarget_matrix = hcat(Xtarget, Xtarget_index)
            Xtarget_matrix = Xtarget_matrix[sortperm(Xtarget_matrix[:, size(Xtarget,2)+1], rev = false), :]
            KZCT_algorithm.Xtarget = Xtarget_matrix[:,1:size(Xtarget,2)]

        elseif targets == :Y
            Ytarget = vcat(KZCT_algorithm.Dˢ_union_E.Ytarget, KZCT_algorithm.D_exluding_Dˢ_union_Ε.Ytarget)
            Ytarget_index = vcat(KZCT_algorithm.Dˢ_union_E.indexDMU, KZCT_algorithm.D_exluding_Dˢ_union_Ε.indexDMU)
            Ytarget_matrix = hcat(Ytarget, Ytarget_index)
            Ytarget_matrix = Ytarget_matrix[sortperm(Ytarget_matrix[:, size(Ytarget,2)+1], rev = false), :]
            KZCT_algorithm.Ytarget = Ytarget_matrix[:,1:size(Ytarget,2)]
        else 
            throw(ArgumentError("`targets` must be :X or :Y"));
        end
    else 
        if targets == :X
            Xtarget = vcat(KZCT_algorithm.Dˢ.Xtarget, KZCT_algorithm.D_exluding_Dˢ.Xtarget)
            Xtarget_index = vcat(KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
            Xtarget_matrix = hcat(Xtarget, Xtarget_index)
            Xtarget_matrix = Xtarget_matrix[sortperm(Xtarget_matrix[:, size(Xtarget,2)+1], rev = false), :]
            KZCT_algorithm.Xtarget = Xtarget_matrix[:,1:size(Xtarget,2)]

        elseif targets == :Y
            Ytarget = vcat(KZCT_algorithm.Dˢ.Ytarget, KZCT_algorithm.D_exluding_Dˢ.Ytarget)
            Ytarget_index = vcat(KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.D_exluding_Dˢ.indexDMU)
            Ytarget_matrix = hcat(Ytarget, Ytarget_index)
            Ytarget_matrix = Ytarget_matrix[sortperm(Ytarget_matrix[:, size(Ytarget,2)+1], rev = false), :]
            KZCT_algorithm.Ytarget = Ytarget_matrix[:,1:size(Ytarget,2)]
        else 
            throw(ArgumentError("`targets` must be :X or :Y"));
        end
    end
end

"""
    getresults!(KZCT_algorithm, Exterior_presence)

    Get the results in the proper order of the sample D 
"""
function getresults!(KZCT_algorithm::KZCTAlgorithm, Exterior_presence::Bool)

    # effi
    getscores!(KZCT_algorithm, Exterior_presence)

    # Lambda
    getlambdas!(KZCT_algorithm, Exterior_presence)

    # slacks 
    getslacks!(KZCT_algorithm, :X, Exterior_presence)
    getslacks!(KZCT_algorithm, :Y, Exterior_presence)

    # targets 
    gettargets!(KZCT_algorithm, :X, Exterior_presence)
    gettargets!(KZCT_algorithm, :Y, Exterior_presence)


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
        namesDMU::Union{Vector{String},Nothing} = nothing,
        optimizer::Union{DEAOptimizer,Nothing} = nothing)
    

    KZCT_algorithm = KZCTAlgorithm(n = size(X, 1), m = size(X,2), s = size(Y,2),
                                    orient = orient, rts = rts, slack = slack, disposX = disposX, disposY = disposY,
                                    dmunames = namesDMU, optimizer = optimizer)

    # Initialize the sample 
    KZCT_algorithm.D = Subset(X = X, Y = Y, indexDMU = vec([i for i in 1:size(X,1)]))

    # Create and fill the subsample
    initialsubset!(KZCT_algorithm)

    # Find the best-practices B_S in D_S
    bestpracticesfinder!(KZCT_algorithm, :Dˢ)

    # Find exterior DMUs in D excluding D^S respect to the hull of B^S 
    exteriorsfinder!(KZCT_algorithm, :Bˢ)

    # If E is not empty, then we have to find best practice DMUs in D^S Union E, otherwise, F = B_S and all DMUs are already evaluated 
    if size(KZCT_algorithm.Ε.X,1)>0

        KZCT_algorithm.Dˢ_union_E = unionsubset(X,Y, KZCT_algorithm.Dˢ.indexDMU, KZCT_algorithm.Ε.indexDMU)
        bestpracticesfinder!(KZCT_algorithm, :Dˢ_union_E)
        KZCT_algorithm.D_exluding_Dˢ_union_Ε = excludingsubset(X,Y, KZCT_algorithm.Dˢ_union_E.indexDMU)
        exteriorsfinder!(KZCT_algorithm, :F) 
        getresults!(KZCT_algorithm, size(KZCT_algorithm.Ε.X,1)>0)
 
    else
        KZCT_algorithm.F = KZCT_algorithm.Bˢ
        getresults!(KZCT_algorithm, size(KZCT_algorithm.Ε.X,1)>0)
    end 
 
    n = KZCT_algorithm.n 
    m = KZCT_algorithm.m 
    s = KZCT_algorithm.s 
    orient = KZCT_algorithm.orient 
    rts = KZCT_algorithm.rts 
    disposX = KZCT_algorithm.disposX
    dispoY = KZCT_algorithm.disposY 
    names = KZCT_algorithm.dmunames
    effi = KZCT_algorithm.eff
    slackX = KZCT_algorithm.slackX 
    slackY = KZCT_algorithm.slackY 
    lambdaeff = KZCT_algorithm.lambda 
    Xtarget = KZCT_algorithm.Xtarget
    Ytarget = KZCT_algorithm.Ytarget

    return RadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end
