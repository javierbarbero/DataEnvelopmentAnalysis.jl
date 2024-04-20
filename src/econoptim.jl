# This file contains functions for economic optimization problems.
# These functions are intended for internal use and package extensions.
function deamincost(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}; rts::Symbol = :VRS, dispos::Symbol = :Strong,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nw != nx
        throw(DimensionMismatch("number of rows in W and X ($nw, $nx) are not equal"));
    end
    if mw != m
        throw(DimensionMismatch("number of columns in W and X ($mw, $m) are not equal"));
    end

    if dispos != :Strong && dispos != :Weak
        throw(ArgumentError("`disposX` must be :Strong or :Weak"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    # Compute efficiency for each DMU
    n = nx

    Xtarget = zeros(n,m)
    Ytarget = Y[:,:]
    clambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        y0 = Y[i,:]
        w0 = W[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, lambda[1:n] >= 0)

        @objective(deamodel, Min, sum(w0[j] .* Xeff[j] for j in 1:m))

        @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) <= Xeff[j])

        if dispos == :Strong
            @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) >= y0[j])
        elseif dispos == :Weak
            @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) == y0[j])
        end

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        elseif rts == :FDH
            @constraint(deamodel, sum(lambda) == 1)
            set_binary.(lambda[1:n])
        else
            throw(ArgumentError("`rts` must be :CRS, :VRS or :FDH"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        Xtarget[i,:]  = JuMP.value.(Xeff)
        clambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    return Xtarget::Matrix, clambdaeff::SparseMatrixCSC{Float64, Int64}

end


function deamaxrevenue(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    P::Union{Matrix,Vector}; rts::Symbol = :VRS, dispos::Symbol = :Strong,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    np, sp = size(P, 1), size(P, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if np != ny
        throw(DimensionMismatch("number of rows in P and Y ($np, $ny) are not equal"));
    end
    if sp != s
        throw(DimensionMismatch("number of columns in P and Y ($sp, $s) are not equal"));
    end

    if dispos != :Strong && dispos != :Weak
        throw(ArgumentError("`disposY` must be :Strong or :Weak"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end    

    # Compute efficiency for each DMU
    n = nx

    Xtarget = X[:,:]
    Ytarget = zeros(n,s)
    rlambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, Yeff[1:s])
        @variable(deamodel, lambda[1:n] >= 0)

        @objective(deamodel, Max, sum(p0[j] .* Yeff[j] for j in 1:s))

        if dispos == :Strong
            @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) <= x0[j])
        elseif dispos == :Weak
            @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) == x0[j])
        end

        @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) >= Yeff[j])

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        elseif rts == :FDH
            @constraint(deamodel, sum(lambda) == 1)
            set_binary.(lambda[1:n])
        else
            throw(ArgumentError("`rts` must be :CRS, :VRS or :FDH"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        Ytarget[i,:]  = JuMP.value.(Yeff)
        rlambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    return Ytarget::Matrix, rlambdaeff::SparseMatrixCSC{Float64, Int64}

end


function deamaxprofit(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}, P::Union{Matrix,Vector}; rts::Symbol = :VRS,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)
    np, sp = size(P, 1), size(P, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nw != nx
        throw(DimensionMismatch("number of rows in W and X ($nw, $nx) are not equal"));
    end
    if np != ny
        throw(DimensionMismatch("number of rows in P and Y ($np, $ny) are not equal"));
    end
    if mw != m
        throw(DimensionMismatch("number of columns in W and X ($mw, $m) are not equal"));
    end
    if sp != s
        throw(DimensionMismatch("number of columns in P and Y ($sp, $s) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end   

    # Compute efficiency for each DMU
    n = nx

    Xtarget = zeros(n,m)
    Ytarget = zeros(n,s)
    plambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, Yeff[1:s])
        @variable(deamodel, lambda[1:n] >= 0)

        @objective(deamodel, Max, (sum(p0[j] .* Yeff[j] for j in 1:s)) - (sum(w0[j] .* Xeff[j] for j in 1:m)))

        @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) <= Xeff[j])
        @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) >= Yeff[j])

        # Add return to scale constraints
        if rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        elseif rts == :FDH
            @constraint(deamodel, sum(lambda) == 1)
            set_binary.(lambda[1:n])
        else
            throw(ArgumentError("`rts` must be :VRS or :FDH"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        Xtarget[i,:]  = JuMP.value.(Xeff)
        Ytarget[i,:]  = JuMP.value.(Yeff)
        plambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    return Xtarget::Matrix, Ytarget::Matrix, plambdaeff::SparseMatrixCSC{Float64, Int64}

end
