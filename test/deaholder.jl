# Tests for Hölder DEA Models
@testset "HolderDEAModel" begin

    # ------------------
    # Input oriented
    # ------------------
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6 8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    # Hölder l = 1
    holderio1 = deaholder(X, Y, l = 1, orient = :Input, rts = :VRS)

    @test typeof(holderio1) == HolderL1DEAModel

    @test nobs(holderio1) == 8
    @test ninputs(holderio1) == 2
    @test noutputs(holderio1) == 1
    @test efficiency(holderio1) ≈ [0; 0; 0; 2.0; 4.0; 0.0; 1.0; 0.6]
    @test (efficiency(holderio1, :min) == [1; 1; 1; 2; 1; 2; 1; 1]) || (efficiency(holderio1, :min) == [1; 2; 1; 2; 1; 2; 1; 1])
    @test convert(Matrix, peers(holderio1)) ≈ 
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0]
    @test slacks(holderio1, :X) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 1.0; 2.0 0.0; 0.0 1.0; 0.0 4.0]
    @test slacks(holderio1, :Y) ≈ zeros(8,1)

    @test efficiency(deaholder(targets(holderio1, :X), targets(holderio1, :Y), l = 1, orient = :Input, rts = :VRS)) ≈ zeros(8,1) atol = 1e-10

    # Test no slacks
    holderio1noslacks = deaholder(X, Y, l = 1, orient = :Input, rts = :VRS, slack = false)
    @test efficiency(holderio1noslacks) == efficiency(holderio1)
    @test isempty(slacks(holderio1noslacks, :X)) == 1
    @test isempty(slacks(holderio1noslacks, :Y)) == 1  
    
    # Weighted (weakly)
    holderio1weight = deaholder(X, Y, l = 1, orient = :Input, rts = :VRS, weight = true)
    @test efficiency(holderio1weight) ≈ [0; 0; 0; 0.625; 0.8; 0.0; 0.5; 0.375]

    # Hölder l = 2
    logs, value = Test.collect_test_logs() do
        holderio2 = deaholder(X, Y, l = 2, orient = :Input, rts = :VRS, optimizer = DEAOptimizer(:NLP))
    end

    @test typeof(holderio2) == HolderL2DEAModel

    @test nobs(holderio2) == 8
    @test ninputs(holderio2) == 2
    @test noutputs(holderio2) == 1

    if (Base.find_package("Gurobi") !== nothing)
        using Gurobi
        @info ("Testing Hölder L2 with Gurobi")

        holderio2 = deaholder(X, Y, l = 2, orient = :Input, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holderio2) ≈ [0.0; 0.0; 0.0; 1.788854; 4.0; 0.0; 1.0; 0.6] atol = 1e-5
        @test slacks(holderio2, :X) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 1.0; 0.0 4.0]
        @test slacks(holderio2, :Y) ≈ zeros(8,1)
        @test convert(Matrix, peers(holderio2)) ≈
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.4  0.0  0.6  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0]

        @test efficiency(deaholder(targets(holderio2, :X), targets(holderio2, :Y), l = 2, orient = :Input, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))) ≈ zeros(8,1)  
        
        # Weighted (weakly)
        holderio2weight = deaholder(X, Y, l = 2, orient = :Input, rts = :VRS, weight = true, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holderio2weight) ≈ [0.0; 0.0; 0.0; 0.554700; 0.8; 0.0; 0.468521; 0.375] atol = 1e-5
    else
        # Function runs without error
        @test occursin("Model solved with Ipopt.Optimizer is innacuate. Use a solver that supports SOS1 constraints.", string(logs))

        # Weighted (weakly)
        logs, value = Test.collect_test_logs() do
            holderio2weight = deaholder(X, Y, l = 2, orient = :Input, rts = :VRS, weight = true, optimizer = DEAOptimizer(:NLP))
        end
    end    

    # Hölder l = Inf
    holderioinf = deaholder(X, Y, l = Inf, orient = :Input, rts = :VRS)

    @test typeof(holderioinf) == HolderLInfDEAModel

    @test nobs(holderioinf) == 8
    @test ninputs(holderioinf) == 2
    @test noutputs(holderioinf) == 1
    @test efficiency(holderioinf) ≈ [0.0; 0.0; 0.0; 4/3; 3.0; 0.0; 1.0; 0.6]
    @test convert(Matrix, peers(holderioinf)) ≈ 
            [ 1.0       0.0  0.0       0.0  0.0  0.0  0.0  0.0
            0.0       1.0  0.0       0.0  0.0  0.0  0.0  0.0
            0.0       0.0  1.0       0.0  0.0  0.0  0.0  0.0
            2/3       0.0  1/3       0.0  0.0  0.0  0.0  0.0
            1.0       0.0  0.0       0.0  0.0  0.0  0.0  0.0
            0.0       0.0  1.0       0.0  0.0  0.0  0.0  0.0
            0.0       1.0  0.0       0.0  0.0  0.0  0.0  0.0
            0.0       1.0  0.0       0.0  0.0  0.0  0.0  0.0]
    @test slacks(holderioinf, :X) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 2.0 0.0; 0.0 0.0; 0.0 3.4]
    @test slacks(holderioinf, :Y) ≈ zeros(8,1)

    @test efficiency(deaholder(targets(holderioinf, :X), targets(holderioinf, :Y), l = 1, orient = :Input, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    holderioinfnoslacks = deaholder(X, Y, l = Inf, orient = :Input, rts = :VRS, slack = false)
    @test efficiency(holderioinfnoslacks) == efficiency(holderioinf)
    @test isempty(slacks(holderioinfnoslacks, :X)) == 1
    @test isempty(slacks(holderioinfnoslacks, :Y)) == 1   

    # Weighted (weakly)
    holderioinfweight = deaholder(X, Y, l = Inf, orient = :Input, rts = :VRS, weight = true)
    @test efficiency(holderioinfweight) ≈ [0; 0; 0; 0.4; 0.6; 0.0; 1/3; 0.375]

    # Print
    show(IOBuffer(), holderio1)
    show(IOBuffer(), holderio1noslacks)
    show(IOBuffer(), holderio1weight)
    show(IOBuffer(), holderio2)
    show(IOBuffer(), holderio2weight)
    show(IOBuffer(), holderioinf)
    show(IOBuffer(), holderioinfnoslacks)
    show(IOBuffer(), holderioinfweight)

    # ------------------
    # Output oriented
    # ------------------
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]    

    # Hölder l = 1
    holderoo1 = deaholder(X, Y, l = 1, orient = :Output, rts = :VRS)

    @test typeof(holderoo1) == HolderL1DEAModel

    @test nobs(holderoo1) == 8
    @test ninputs(holderoo1) == 1
    @test noutputs(holderoo1) == 2
    @test efficiency(holderoo1) ≈ [0; 0; 0; 3.0; 5.0; 0.0; 2.0; 3.0]
    @test (efficiency(holderoo1, :min) == [2; 2; 2; 3; 2; 2; 2; 3]) || (efficiency(holderoo1, :min) ==  [3; 2; 3; 3; 2; 2; 2; 3])
    @test convert(Matrix, peers(holderoo1)) ≈ 
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0]
    @test slacks(holderoo1, :X) ≈ zeros(8,1)
    @test slacks(holderoo1, :Y) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 2.0; 0.0 0.0; 2.5 0.0]

    @test efficiency(deaholder(targets(holderoo1, :X), targets(holderoo1, :Y), l = 1, orient = :Output, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    holderoo1noslacks = deaholder(X, Y, l = 1, orient = :Output, rts = :VRS, slack = false)
    @test efficiency(holderoo1noslacks) == efficiency(holderoo1)
    @test isempty(slacks(holderoo1noslacks, :X)) == 1
    @test isempty(slacks(holderoo1noslacks, :Y)) == 1  
    
    # Weighted (weakly)
    holderoo1weight = deaholder(X, Y, l = 1, orient = :Output, rts = :VRS, weight = true)
    @test efficiency(holderoo1weight) ≈ [0; 0; 0; 0.6; 5/3; 0.0; 1/3; 0.6]

    # Hölder l = 2
    logs, value = Test.collect_test_logs() do
        holderoo2 = deaholder(X, Y, l = 2, orient = :Output, rts = :VRS, optimizer = DEAOptimizer(:NLP))
    end

    @test typeof(holderoo2) == HolderL2DEAModel

    @test nobs(holderoo2) == 8
    @test ninputs(holderoo2) == 1
    @test noutputs(holderoo2) == 2

    if (Base.find_package("Gurobi") !== nothing)
        using Gurobi
        @info ("Testing Hölder L2 with Gurobi")

        holderoo2 = deaholder(X, Y, l = 2, orient = :Output, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holderoo2) ≈ [0.0; 0.0; 0.0; 3.0; 5.0; 0.0; 1.897367; 3.0] atol = 1e-5
        @test slacks(holderoo2, :X) ≈ zeros(8,1)
        @test slacks(holderoo2, :Y) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 1.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 0.0; 2.5 0.0]
        @test convert(Matrix, peers(holderoo2)) ≈
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
            0.2  0.0  0.8  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0]

        @test efficiency(deaholder(targets(holderoo2, :X), targets(holderoo2, :Y), l = 2, orient = :Output, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))) ≈ zeros(8,1) atol = 1e-5 
        
        # Weighted (weakly)
        holderoo2weight = deaholder(X, Y, l = 2, orient = :Output, rts = :VRS, weight = true, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holderoo2weight) ≈ [0.0; 0.0; 0.0; 0.6; 5/3; 0.0; 0.325396; 0.6] atol = 1e-5
    else
        # Function runs without error
        @test occursin("Model solved with Ipopt.Optimizer is innacuate. Use a solver that supports SOS1 constraints.", string(logs))

        # Weighted (weakly)
        logs, value = Test.collect_test_logs() do
            holderoo2weight = deaholder(X, Y, l = 2, orient = :Output, rts = :VRS, weight = true, optimizer = DEAOptimizer(:NLP))
        end
    end    

    # Hölder l = Inf
    holderooinf = deaholder(X, Y, l = Inf, orient = :Output, rts = :VRS)

    @test typeof(holderooinf) == HolderLInfDEAModel

    @test nobs(holderooinf) == 8
    @test ninputs(holderooinf) == 1
    @test noutputs(holderooinf) == 2
    @test efficiency(holderooinf) ≈ [0.0; 0.0; 0.0; 2.5; 4.0; 0.0; 1.5; 2.875]
    @test convert(Matrix, peers(holderooinf)) ≈ 
            [ 1.0    0.0    0.0  0.0  0.0  0.0  0.0  0.0
            0.0    1.0    0.0  0.0  0.0  0.0  0.0  0.0
            0.0    0.0    1.0  0.0  0.0  0.0  0.0  0.0
            0.5    0.5    0.0  0.0  0.0  0.0  0.0  0.0
            1.0    0.0    0.0  0.0  0.0  0.0  0.0  0.0
            0.0    0.0    1.0  0.0  0.0  0.0  0.0  0.0
            0.5    0.0    0.5  0.0  0.0  0.0  0.0  0.0
            0.125  0.875  0.0  0.0  0.0  0.0  0.0  0.0]
    @test slacks(holderooinf, :X) ≈ zeros(8,1)
    @test slacks(holderooinf, :Y) ≈ [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 2.0; 0.0 0.0; 0.0 0.0]

    @test efficiency(deaholder(targets(holderooinf, :X), targets(holderooinf, :Y), l = 1, orient = :Output, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    holderooinfnoslacks = deaholder(X, Y, l = Inf, orient = :Output, rts = :VRS, slack = false)
    @test efficiency(holderooinfnoslacks) == efficiency(holderooinf)
    @test isempty(slacks(holderooinfnoslacks, :X)) == 1
    @test isempty(slacks(holderooinfnoslacks, :Y)) == 1   

    # Weighted (weakly)
    holderooinfweight = deaholder(X, Y, l = Inf, orient = :Output, rts = :VRS, weight = true)
    @test efficiency(holderooinfweight) ≈ [0.0; 0.0; 0.0; 5/9; 4/3; 0.0; 3/11; 0.6]

    # Print
    show(IOBuffer(), holderoo1)
    show(IOBuffer(), holderoo1noslacks)
    show(IOBuffer(), holderoo1weight)
    show(IOBuffer(), holderoo2)
    show(IOBuffer(), holderoo2weight)
    show(IOBuffer(), holderooinf)
    show(IOBuffer(), holderooinfnoslacks)
    show(IOBuffer(), holderooinfweight)

    # ------------------
    # Graph
    # ------------------
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]    

    # Hölder l = 1
    holdergr1 = deaholder(X, Y, l = 1, orient = :Graph, rts = :VRS)

    @test typeof(holdergr1) == HolderL1DEAModel

    @test nobs(holdergr1) == 8
    @test ninputs(holdergr1) == 1
    @test noutputs(holdergr1) == 1
    @test efficiency(holdergr1) ≈ [0; 0; 0; 0.0; 3.0; 2.0; 0.0; 6.0]
    @test (efficiency(holdergr1, :min) == [1; 1; 1; 1; 1; 2; 2; 2]) || (efficiency(holdergr1, :min)  == [2; 1; 1; 2; 1; 2; 2; 2])
    @test convert(Matrix, peers(holdergr1)) ≈ 
            [1.0  0.0  0.0    0.0    0.0  0.0  0.0  0.0
            0.0  1.0  0.0    0.0    0.0  0.0  0.0  0.0
            0.0  0.0  1.0    0.0    0.0  0.0  0.0  0.0
            0.0  0.0  0.0    1.0    0.0  0.0  0.0  0.0
            0.5  0.5  0.0    0.0    0.0  0.0  0.0  0.0
            0.0  0.0  0.0    1.0    0.0  0.0  0.0  0.0
            0.0  0.0  0.0    1.0    0.0  0.0  0.0  0.0
            0.0  0.0  0.647  0.353  0.0  0.0  0.0  0.0]
    @test slacks(holdergr1, :X) ≈ [0; 0; 0; 0.0; 0.0; 2.0; 2.0; 0.0]
    @test slacks(holdergr1, :Y) ≈ zeros(8,1) atol = 1e-5

    @test efficiency(deaholder(targets(holdergr1, :X), targets(holdergr1, :Y), l = 1, orient = :Graph, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    holdergr1noslacks = deaholder(X, Y, l = 1, orient = :Graph, rts = :VRS, slack = false)
    @test efficiency(holdergr1noslacks) == efficiency(holdergr1)
    @test isempty(slacks(holdergr1noslacks, :X)) == 1
    @test isempty(slacks(holdergr1noslacks, :Y)) == 1  
        
    # Weighted (weakly)
    holdergr1weight = deaholder(X, Y, l = 1, orient = :Graph, rts = :VRS, weight = true)
    @test efficiency(holdergr1weight) ≈ [0; 0; 0; 0.0; 0.5; 2/7; 0.0; 0.715629] atol = 1e-5
    
    # Hölder l = 2
    logs, value = Test.collect_test_logs() do
        holdergr2 = deaholder(X, Y, l = 2, orient = :Graph, rts = :VRS, optimizer = DEAOptimizer(:NLP))
    end

    @test typeof(holdergr2) == HolderL2DEAModel

    @test nobs(holdergr2) == 8
    @test ninputs(holdergr2) == 1
    @test noutputs(holdergr2) == 1
    
    if (Base.find_package("Gurobi") !== nothing)
        using Gurobi
        @info ("Testing Hölder L2 with Gurobi")

        holdergr2 = deaholder(X, Y, l = 2, orient = :Graph, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holdergr2) ≈ [0.0; 0.0; 0.0; 0.0; 2.683282; 2.0; 0.0; 5.3648] atol = 1e-5
        @test slacks(holdergr2, :X) ≈ [0.0; 0.0; 0.0; 0.0; 0.0; 2.0; 0.0; 0.0]
        @test slacks(holdergr2, :Y) ≈ [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] 
        @test convert(Matrix, peers(holdergr2)) ≈
            [1.0  0.0      0.0      0.0  0.0  0.0  0.0  0.0
            0.0  1.0      0.0      0.0  0.0  0.0  0.0  0.0
            0.0  0.0      1.0      0.0  0.0  0.0  0.0  0.0
            0.0  0.0      0.0      1.0  0.0  0.0  0.0  0.0
            0.2  0.8      0.0      0.0  0.0  0.0  0.0  0.0
            0.0  0.0      0.0      1.0  0.0  0.0  0.0  0.0
            0.0  0.0      0.0      0.0  0.0  0.0  1.0  0.0
            0.0  0.45172  0.54828  0.0  0.0  0.0  0.0  0.0]

        @test efficiency(deaholder(targets(holdergr2, :X), targets(holdergr2, :Y), l = 2, orient = :Graph, rts = :VRS, optimizer = DEAOptimizer(Gurobi.Optimizer))) ≈ zeros(8,1) atol = 1e-5 
        
        # Weighted (weakly)
        holdergr2weight = deaholder(X, Y, l = 2, orient = :Graph, rts = :VRS, weight = true, optimizer = DEAOptimizer(Gurobi.Optimizer))
        @test efficiency(holdergr2weight) ≈  [0.0; 0.0; 0.0; 0.0; 0.485071; 2/7; 0.0; 0.710103] atol = 1e-5
    else
        # Function runs without error
        @test occursin("Model solved with Ipopt.Optimizer is innacuate. Use a solver that supports SOS1 constraints.", string(logs))

        # Weighted (weakly)
        logs, value = Test.collect_test_logs() do
            holdergr2weight = deaholder(X, Y, l = 2, orient = :Graph, rts = :VRS, weight = true, optimizer = DEAOptimizer(:NLP))
        end
    end    

    # Hölder l = Inf
    holdergrinf = deaholder(X, Y, l = Inf, orient = :Graph, rts = :VRS)

    @test typeof(holdergrinf) == HolderLInfDEAModel

    @test nobs(holdergrinf) == 8
    @test ninputs(holdergrinf) == 1
    @test noutputs(holdergrinf) == 1
    @test efficiency(holdergrinf) ≈ [0.0; 0.0; 0.0; 0.0; 2.0; 2.0; 0.0; 3.832]
    @test convert(Matrix, peers(holdergrinf)) ≈ 
            [ 1.0  0.0    0.0    0.0  0.0  0.0  0.0  0.0
            0.0  1.0    0.0    0.0  0.0  0.0  0.0  0.0
            0.0  0.0    1.0    0.0  0.0  0.0  0.0  0.0
            0.0  0.0    0.0    1.0  0.0  0.0  0.0  0.0
            0.0  1.0    0.0    0.0  0.0  0.0  0.0  0.0
            0.0  0.0    0.0    1.0  0.0  0.0  0.0  0.0
            0.0  0.0    0.0    1.0  0.0  0.0  0.0  0.0
            0.0  0.605  0.395  0.0  0.0  0.0  0.0  0.0]
    @test slacks(holdergrinf, :X) ≈ [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 2.0; 0.0] atol = 1e-5
    @test slacks(holdergrinf, :Y) ≈ [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] atol = 1e-5

    @test efficiency(deaholder(targets(holdergrinf, :X), targets(holdergrinf, :Y), l = 1, orient = :Graph, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    holdergrinfnoslacks = deaholder(X, Y, l = Inf, orient = :Graph, rts = :VRS, slack = false)
    @test efficiency(holdergrinfnoslacks) == efficiency(holdergrinf)
    @test isempty(slacks(holdergrinfnoslacks, :X)) == 1
    @test isempty(slacks(holdergrinfnoslacks, :Y)) == 1   

    # Weighted (weakly)
    holdergrinfweight = deaholder(X, Y, l = Inf, orient = :Graph, rts = :VRS, weight = true)
    @test efficiency(holdergrinfweight) ≈ [0.0; 0.0; 0.0; 0.0; 0.4; 5/21; 0.0; 0.636115] atol = 1e-5

    # Print
    show(IOBuffer(), holdergr1)
    show(IOBuffer(), holdergr1noslacks)
    show(IOBuffer(), holdergr1weight)
    show(IOBuffer(), holdergr2)
    show(IOBuffer(), holdergr2weight)
    show(IOBuffer(), holdergrinf)
    show(IOBuffer(), holdergrinfnoslacks)
    show(IOBuffer(), holdergrinfweight)

    # ------------------
    # Test errors
    # ------------------
    @test_throws DimensionMismatch deaholder([1; 2 ; 3], [4 ; 5], l = 1) #  Different number of observations
    @test_throws DimensionMismatch deaholder([1; 2], [4 ; 5], l = 1, Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deaholder([1 1; 2 2], [4 4; 5 5], l = 1, Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deaholder([1 1; 2 2], [4 4; 5 5], l = 1, Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deaholder([1; 2; 3], [4; 5; 6], l = 1, orient = :Error) # Invalid orientation
    @test_throws ArgumentError deaholder([1; 2; 3], [4; 5; 6], l = 1, rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError deaholder([1; 2; 3], [4; 5; 6], l = 2, rts = :Error) # CRS in L2
    @test_throws ArgumentError deaholder([1; 2; 3], [4; 5; 6], l = 5) # Invalid r
    @test_throws ArgumentError deaholder([1; 2; 3], [4; 5; 6], l = 2, rts = :VRS) # No solver
    @test_throws ArgumentError efficiency(holderio1, :Error) # Invalid efficiency type

end
