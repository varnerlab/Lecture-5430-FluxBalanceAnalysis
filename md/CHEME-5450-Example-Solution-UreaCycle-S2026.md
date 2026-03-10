# Problem Set 2 (PS2): Flux Balance Analysis of the Urea Cycle in HL-60 Cells
In problem set 2 (PS2), we will explore the urea cycle in HL-60 cells using flux balance analysis. The [urea cycle](https://www.kegg.jp/pathway/hsa00220) is a crucial metabolic pathway that converts toxic ammonia into urea for excretion. While the urea cycle's role in [HL-60 cells, a human promyelocytic leukemia cell line](https://www.atcc.org/products/ccl-240?matchtype=b&network=g&device=c&adposition=&keyword=hl60%20cell%20line%20atcc&gad_source=1&gbraid=0AAAAADR6fpoOXsp8U8fXLd_E6sLTcwv24&gclid=CjwKCAiA5eC9BhAuEiwA3CKwQm0C1oE5_JjTpJ24VnTjZUZQVLivpPxmufDo7HdH5v3hN1XKnEf3ExoCvhwQAvD_BwE), is not directly established, these cells exhibit alterations in protein levels and proliferation rates when exposed to various compounds, which may indirectly affect nitrogen metabolism and related pathways.

### Tasks
In PS2, we'll construct [a simplified model of the urea cycle](https://github.com/varnerlab/CHEME-5450-Lectures-Spring-2025/blob/main/lectures/week-5/L5c/docs/figs/Fig-Urea-cycle-Schematic.pdf), determine the reversibility of the reactions, estimates for the flux bounds, and then compute the optimal flux distribution that maximizes Urea production. Start with the setup section, and work your way through the notebook. `TODO` statements/comments indicate that you need to do something.

### References
1. [Al-Otaibi NAS, Cassoli JS, Martins-de-Souza D, Slater NKH, Rahmoune H. Human leukemia cells (HL-60) proteomic and biological signatures underpinning cryo-damage are differentially modulated by novel cryo-additives. Gigascience. 2019 Mar 1;8(3):giy155. doi: 10.1093/gigascience/giy155. PMID: 30535373; PMCID: PMC6394207.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6394207/)
2. [Figarola JL, Weng Y, Lincoln C, Horne D, Rahbar S. Novel dichlorophenyl urea compounds inhibit proliferation of human leukemia HL-60 cells by inducing cell cycle arrest, differentiation and apoptosis. Invest New Drugs. 2012 Aug;30(4):1413-25. doi: 10.1007/s10637-011-9711-8. Epub 2011 Jul 5. PMID: 21728022.](https://pubmed.ncbi.nlm.nih.gov/21728022/)
3. [Caldwell RW, Rodriguez PC, Toque HA, Narayanan SP, Caldwell RB. Arginase: A Multifaceted Enzyme Important in Health and Disease. Physiol Rev. 2018 Apr 1;98(2):641-665. doi: 10.1152/physrev.00037.2016. PMID: 29412048; PMCID: PMC5966718.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5966718/)

## Setup, Data, and Prerequisites
We set up the computational environment by including the `Include.jl` file and loading any needed resources.

> The [`include(...)` command](https://docs.julialang.org/en/v1/base/base/#include) evaluates the contents of the input source file, `Include.jl`, in the notebook's global scope. The `Include.jl` file sets paths, loads required external packages, custom types, and functions used in this exercise. It checks for a `Manifest.toml` file; if one is found, packages are loaded from it. Otherwise, packages are downloaded and loaded.

Let's set up our code environment:


```julia
include("Include.jl");
```

__Build the model__. To store all the problem data, we created [the `MyPrimalFluxBalanceAnalysisCalculationModel` type](src/Types.jl). Let's build one of these objects for our problem and store it in the `model::MyPrimalFluxBalanceAnalysisCalculationModel` variable. We also return the `rd::Dict{String, String}` dictionary, which maps the reaction name field (key) to the reaction string (value).
* __Builder (or factory) pattern__: For all custom types that we make, we'll use something like [the builder software pattern](https://en.wikipedia.org/wiki/Builder_pattern) to construct and initialize these objects. The calling syntax will be the same for all types: [a `build(...)` method](src/Factory.jl) will take the kind of thing we want to build in the first argument, and the data needed to make that type as [a `NamedTuple` instance](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) in the second argument.
* __What's the story with the `let` block__? A [let block](https://docs.julialang.org/en/v1/manual/variables-and-scoping/#Let-Blocks) creates a new hard scope and new variable bindings each time they run. Thus, they act like a private scratch space, where data comes in (is captured by the block), but only what we want to be exposed comes out. 


```julia
model, rd = let

    # first, load the reaction file - and process it
    listofreactions = read_reaction_file(joinpath(_PATH_TO_DATA, "Network.net")); # load the reactions from the VFF reaction file
    S, species, reactions, rd = build_stoichiometric_matrix(listofreactions); # Builds the stochiometric matrix, species list, and the reactions list
    boundsarray = build_default_bounds_array(listofreactions); # Builds a default bounds model using the flat file flags

    # build the FBA model -
    model = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
        S = S, # stoichiometric matrix
        fluxbounds = boundsarray, # these are the *default* bounds, we'll need to update with new info if we have it
        species = species, # list of species. The rows of S are in this order
        reactions = reactions, # list of reactions. The cols of S are in this order
        objective = length(reactions) |> R -> zeros(R), # this is empty, we'll need to set this
    ));

    # return -
    model, rd
end;
```

## Update the flux bounds
The flux bounds are important constraints in flux balance analysis calculations and the convex decomposition of the stoichiometric array. Beyond their role in the flux estimation problem, the flux bounds are _integrative_, i.e., these constraints integrate many types of genetic and biochemical information into the problem. A general model for these bounds is given by:
$$
\begin{align*}
-\delta_{j}\underbrace{\left[{V_{max,j}^{\circ}}\left(\frac{e}{e^{\circ}}\right)\theta_{j}\left(\dots\right){f_{j}\left(\dots\right)}\right]}_{\text{reverse: other functions or parameters?}}\leq\hat{v}_{j}\leq{V_{max,j}^{\circ}}\left(\frac{e}{e^{\circ}}\right)\theta_{j}\left(\dots\right){f_{j}\left(\dots\right)}
\end{align*}
$$
where $V_{max,j}^{\circ}$ denotes the maximum reaction velocity (units: `flux`) computed at some _characteristic enzyme abundance_. Thus, the maximum reaction velocity is given by:
$$
V_{max,j}^{\circ} \equiv k_{cat,j}^{\circ}e^{\circ}
$$
where $k_{cat,j}$ is the catalytic constant or turnover number for the enzyme (units: `1/time`) and $e^{\circ}$ is a characteristic enzyme abundance (units: `concentration`). The term $\left(e/e^{\circ}\right)$ is a correction to account for the _actual_ enzyme abundance catalyzing the reaction (units: `dimensionless`). The $\theta_{j}\left(\dots\right)\in\left[0,1\right]$ is the current fraction of maximal enzyme activity of enzyme $e$ in reaction $j$. The activity model $\theta_{j}\left(\dots\right)$ describes [allosteric effects](https://en.wikipedia.org/wiki/Allosteric_regulation) on the reaction rate, and is a function of the regulatory and the chemical state of the system, the concentration of substrates, products, and cofactors (units: `dimensionless`).
Finally, the $f_{j}\left(\dots\right)$ is a function describing the substrate (reactants) dependence of the reaction rate $j$ (units: `dimensionless`). 

### PS2
In problem set 2, we use a simplified bounds model to estimate the flux bounds.
We assume that $(e/e^{\circ})\sim{1}$, there are no allosteric inputs $\theta_{j}\left(\dots\right)\sim{1}$, and the substrates are saturating $f_{j}\left(\dots\right)\sim{1}$. 
Then, the flux bounds are given by:
$$
\begin{align*}
-\delta_{j}V_{max,j}^{\circ}\leq{\hat{v}_{j}}\leq{V_{max,j}^{\circ}}
\end{align*}
$$
It is easy to see that the flux bounds are a function of the maximum reaction velocity, the catalytic constant or turnover number, and our assumed value of a characteristic enzyme abundance.

### Estimate the reversibility parameters $\delta_{j}$ for reaction $j$
First, let's estimate the reversibility parameter $\delta_{j}$ for each of the `5` enzyme-catalyzed reactions in our model of the Urea cycle using [eQuilibrator](https://equilibrator.weizmann.ac.il).
* [Beber ME, Gollub MG, Mozaffari D, Shebek KM, Flamholz AI, Milo R, Noor E. eQuilibrator 3.0: a database solution for thermodynamic constant estimation. Nucleic Acids Res. 2022 Jan 7;50(D1): D603-D609. doi: 10.1093/nar/gkab1106. PMID: 34850162; PMCID: PMC8728285.](https://pubmed.ncbi.nlm.nih.gov/34850162/)

Store the results of this analysis in the `reversibility_parameter_dictionary::Dict{String, Bool}` dictionary, which maps the reaction name field (key) to the estimated reversibility parameter (value). Use a $\Delta\bar{{G}}$ threshold value (hyperparameter) of `-10.0 kJ/mol` to determine the reversibility of the reactions. 

`TODO`: Fill in the missing elements in the code block below to complete this task; when using [eQuilibrator](https://equilibrator.weizmann.ac.il), use the EC prefix, i.e., `EC 6.3.4.5` and superscript `m` value:


```julia
reversibility_parameter_dictionary = let

    # initialize -
    ΔḠ = -10.0; # threshold value, units: kJ/mol -
    names = model.reactions; # get an array of the names of the reactions (includes exchange)
    reversibility_parameter_dictionary = Dict{String, Int}();

    # TODO: build a ΔG array for the reactions in the model
    ΔG = [
        -4.3 ; # 1 v₁ EC 6.3.4.5 value: -4.3 ± 2.9 kJ/mol
        -5.5 ; # 2 v₂ EC 4.3.2.1 value: -5.5 ± 5.7 kJ/mol
        -51.0 ; # 3 v₃ EC 3.5.3.1 value: -51 ± 12.4 kJ/mol
        -30.3 ; # 4 v₄ EC 2.1.3.3 value: -30.3 ± 5.7 kJ/mol 
        -1220.2 ; # 5 v₅ EC 1.15.13.39 value: -1220.2 ± 29.6 kJ/mol
    ];

    # compute loop -
    for i ∈ eachindex(names)
        name = names[i]; # get the reaction name for reaction i -

        # check: do we have an exchange flux? If so: skip
        if (contains(name, "b") == true)
            continue;
        end
        
        ΔGᵢ = ΔG[i]; # get the ΔG value for reaction i -
        
        # if ΔGᵢ > ΔḠ (less negative than threshold) → reversible (δ=1)
        # if ΔGᵢ < ΔḠ (more negative than threshold) → irreversible (δ=0)
        δᵢ = sign(ΔGᵢ - ΔḠ) == 1 ? 1 : 0

        # store -
        reversibility_parameter_dictionary[name] = δᵢ; # stores the reversibility parameter with key name
    end

    
    # return -
    reversibility_parameter_dictionary;
end;
```

### Estimate the maximum reaction velocity $V_{max,j}^{\circ}$ for reaction $j$
Next, we'll estimate the maximum reaction velocity $V_{max,j}^{\circ}$ for each of the `5` enzyme-catalyzed reactions in our model of the Urea cycle using [BRENDA](https://www.brenda-enzymes.org/):
* [Antje Chang et al., BRENDA, the ELIXIR core data resource in 2021: new developments and updates, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D498–D508, https://doi.org/10.1093/nar/gkaa1025](https://academic.oup.com/nar/article/49/D1/D498/5992283)

Assume the characteristic enzyme abundance $e^{\circ}\simeq$ 0.01 `μmol/gDW` for all reactions. For missing turnover numbers, use a characteristic value of $k_{cat,j}^{\circ}\simeq$ 10 `1/s`. Store the results of this analysis in the `maximum_reaction_velocity_dictionary::Dict{String, Float64}` dictionary, which maps the reaction name field (key) to the estimated maximum reaction velocity (value). 

`TODO`: Fill in the missing elements in the code block below to complete this task:


```julia
maximum_reaction_velocity_dictionary = let

    # initialize -
    eₒ = 0.01; # characteristic enzyme abundance units: mmol/gDW
    kₒ = 10.0; # characteristic turnover rate units: 1/s (use this if we don't have a specific value from BRENDA)
    names = model.reactions; # get an array of the names of the reactions (includes exchange)
    maximum_reaction_velocity_dictionary = Dict{String, Float64}();

    # TODO: **Update** the kcat array for the reactions in the model with values from BRENDA
    kcat = [
        kₒ ;   # 1 v₁ EC 6.3.4.5 No value, use default
        3.28 ; # 2 v₂ EC 4.3.2.1 value: 3.28 1/s
        190.0 ;# 3 v₃ EC 3.5.3.1 value: 190 1/s
        410.0 ;# 4 v₄ EC 2.1.3.3 value: 410 1/s E. coli value. Good choice?
        kₒ ;   # 5 v₅ EC 1.15.13.39 No value, use default
    ];

    # compute loop -
    for i ∈ eachindex(names)
        name = names[i]; # get the reaction name for reaction i -

        # check: do we have an exchange flux? If so: skip
        if (contains(name, "b") == true)
            continue;
        end
        
        # compute the VMaxᵢ -
        VMaxᵢ = kcat[i]*eₒ;

        # store -
        maximum_reaction_velocity_dictionary[name] = VMaxᵢ; # stores the VMax with key name
    end
    
    # return -
    maximum_reaction_velocity_dictionary;
end;
```

### Update the flux bounds array
Now that we have the reversibility parameters and the maximum reaction velocities, we can update the flux bounds for each reaction in our model of the Urea cycle. 

`TODO`: Fill in the missing elements in the code block below to complete this task:


```julia
fluxbounds = let
    
    fluxbounds = model.fluxbounds;
    names = model.reactions;
    for i ∈ eachindex(names)
        name = names[i]; # get the reaction name for reaction i -
    
        # check: do we have an exchange flux? If so: skip
        if (contains(name, "b") == true)
            continue;
        end
        
        VMax = maximum_reaction_velocity_dictionary[name]; # what is the maximum velocity for this reaction?
        δᵢ = reversibility_parameter_dictionary[name]; # what is the reversibility parameter for this reaction?
    
        # update the bounds: lower bound is -δᵢ*VMax (negative for reversible reactions)
        fluxbounds[i,1] = -δᵢ*VMax; # lower bound
        fluxbounds[i,2] = VMax; # upper bound
    end

    # return -
    fluxbounds;
end;
```


```julia
model.fluxbounds = fluxbounds;
```

### Update the objective function
Finally, let's update the objective function of our model. By default, all the elements of the objective function are set to `0.0`. Update the objective function to _maximize_ the export of  `Urea` from the system. 

`TODO`: Fill in the missing elements in the code block below to complete this task:


```julia
objective = model.objective;
reaction_to_maximize = "b4"; # TODO: specify the reaction to maximize (use the reaction name)
findfirst(x-> x==reaction_to_maximize, model.reactions) |> i -> objective[i] = -1; # why negative 1? because we're maximizing the flux through the reaction
```

## Compute the optimal flux distribution
Let's compute the optimal metabolic distribution $\left\{\hat{v}_{i} \mid i = 1,2,\dots,\mathcal{R}\right\}$ by solving the [linear programming problem](). We solve the optimization problem by passing the `model::MyPrimalFluxBalanceAnalysisCalculationModel` to [the `solve(...)` method](src/Compute.jl). This method returns a `solution::Dict{String, Any}` dictionary containing information about the solution.
* __Why the [try-catch environment](https://docs.julialang.org/en/v1/base/base/#try)__? The [solve(...) method](src/Compute.jl) has an [@assert statement](https://docs.julialang.org/en/v1/base/base/#Base.@assert) to check if the calculation has converged. Thus, the solve method will [throw](https://docs.julialang.org/en/v1/base/base/#Core.throw) an [AssertionError](https://docs.julialang.org/en/v1/base/base/#Core.AssertionError) if the optimization problem fails to converge. To gracefully handle this case, we use a [try-catch construct](https://docs.julialang.org/en/v1/base/base/#try). See the [is_solved_and_feasible method from the JuMP package](https://jump.dev/JuMP.jl/stable/api/JuMP/#JuMP.is_solved_and_feasible) for more information.


```julia
solution = let
    
    solution = nothing; # initialize nothing for the solution
    try
        solution = solve(model); # call the solve method with our problem model -
    catch error
        println("error: $(error)"); # Oooooops! Looks like we have a *major malfunction*, problem didn't solve
    end

    # return solution
    solution
end;
```

### Flux table
`Unhide` the code block below to see how we constructed the flux table using [the `pretty_tables(...)` method exported by the `PrettyTables.jl` package](https://github.com/ronisbr/PrettyTables.jl).
* __Summary__: Each row of the flux table holds information about a reaction in the model. The first column has the reaction name, the second column has the estimated optimal flux value (solution of the FBA problem), the third and fourth columns hold the lower (LB) and upper (UB) for the estimated flux, and the last column has the reaction string for the flux.


```julia
flux_table = let

    # setup -
    S = model.S;
    flux_bounds_array = model.fluxbounds;
    number_of_reactions = size(S, 2);
    flux = solution["argmax"];

    # build DataFrame -
    df = DataFrame(
        reaction  = model.reactions,
        flux      = flux,
        LB        = flux_bounds_array[:, 1],
        UB        = flux_bounds_array[:, 2],
        equation  = [rd[r] for r in model.reactions]
    );

    # let's make a pretty table -
    pretty_table(
        df;
        fit_table_in_display_vertically = false,
        fit_table_in_display_horizontally = false,
        backend = :text,
        alignment = [:l, :r, :r, :r, :l],  # last column left-justified
        table_format = TextTableFormat(borders = text_table_borders__compact)
    );
end
```

     ---------- --------- --------- --------- ----------------------------------------------------------------------------------------------------------------
     [1m reaction [0m [1m    flux [0m [1m      LB [0m [1m      UB [0m [1m equation                                                                                                       [0m
     [90m String   [0m [90m Float64 [0m [90m Float64 [0m [90m Float64 [0m [90m String                                                                                                         [0m
     ---------- --------- --------- --------- ----------------------------------------------------------------------------------------------------------------
      v1          0.0328      -0.1       0.1   M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c = M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
      v2          0.0328   -0.0328    0.0328   M_N-(L-Arginino)succinate_c = M_Fumarate_c+M_L-Arginine_c
      v3          0.0328       0.0       1.9   M_L-Arginine_c+M_H2O_c = M_L-Ornithine_c+M_Urea_c
      v4          0.0328       0.0       4.1   M_Carbamoyl_phosphate_c+M_L-Ornithine_c = M_Orthophosphate_c+M_L-Citrulline_c
      v5             0.0       0.0       0.1   2*M_L-Arginine_c+4*M_Oxygen_c+3*M_NADPH_c+3*M_H_c = 2*M_Nitric_oxide_c+2*M_L-Citrulline_c+3*M_NADP_c+4*M_H2O_c
      b1          0.0328   -1000.0    1000.0   [] = M_Carbamoyl_phosphate_c
      b2          0.0328   -1000.0    1000.0   [] = M_L-Aspartate_c
      b3         -0.0328   -1000.0    1000.0   [] = M_Fumarate_c
      b4         -0.0328   -1000.0    1000.0   [] = M_Urea_c
      b5          0.0328   -1000.0    1000.0   [] = M_ATP_c
      b6         -0.0328   -1000.0    1000.0   [] = M_AMP_c
      b7         -0.0328   -1000.0    1000.0   [] = M_Diphosphate_c
      b8         -0.0328   -1000.0    1000.0   [] = M_Orthophosphate_c
      b9             0.0   -1000.0    1000.0   [] = M_Oxygen_c
      b10            0.0   -1000.0    1000.0   [] = M_NADPH_c
      b11            0.0   -1000.0    1000.0   [] = M_H_c
      b12            0.0   -1000.0    1000.0   [] = M_Nitric_oxide_c
      b13            0.0   -1000.0    1000.0   [] = M_NADP_c
      b14         0.0328   -1000.0    1000.0   [] = M_H2O_c
     ---------- --------- --------- --------- ----------------------------------------------------------------------------------------------------------------



```julia
do_I_see_the_flux_table = true; # TODO: update this flag to {true | false} if the flux table is visible
```

## Discussion
Use you code and simulation results to answer the following questions.

__DQ1__: What is the maximum rate of Urea export from the system using your updated model parameters, and what species are exported or imported into the system to support this production level?


```julia
# Put your answer to DQ1 (either as a commented code cell, or as a markdown cell)
```


```julia
did_I_answer_DQ1 = true; # update to true if answered DQ1 {true | false}
```

__DQ2__: Given your updated model parameters, is there a rate-limiting step controlling the rate of Urea production?


```julia
# Put your answer to DQ2 (either as a commented code cell, or as a markdown cell)
```


```julia
did_I_answer_DQ2 = true; # update to true if answered DQ2 {true | false}
```

__DQ3__: Hypothetically, suppose the teaching team measured the rate of oxygen consumption for the Urea cycle in isolation and found this number to be non-zero. What does that say about the reactions occurring inside the Urea cycle system?


```julia
# Put your answer to DQ3 (either as a commented code cell, or as a markdown cell)
```


```julia
did_I_answer_DQ3 = true; # update to true if answered DQ3 {true | false}
```

## Tests
`Unhide` the code block below (if you are curious) about how we implemented the tests and what we are testing. In these tests, we check values in your notebook and give feedback on which items are correct, missing etc.


```julia
let
    @testset verbose = true "CHEME 5450 problem set 2 test suite" begin
        
        @testset "Setup" begin
            @test isnothing(model) == false
            @test isnothing(rd) == false
            @test isnothing(reversibility_parameter_dictionary) == false
            @test isnothing(maximum_reaction_velocity_dictionary) == false
            @test isnothing(solution) == false

            @test isempty(reversibility_parameter_dictionary) == false
            @test isempty(maximum_reaction_velocity_dictionary) == false

        end

        @testset "Calculation" begin
            @test isempty(solution) == false
            @test do_I_see_the_flux_table == true
        end
        
       @testset "Discussion questions" begin
            @test did_I_answer_DQ1 == true
            @test did_I_answer_DQ2 == true
            @test did_I_answer_DQ3 == true
        end
    end
end;
```

    [0m[1mTest Summary:                       | [22m[32m[1mPass  [22m[39m[36m[1mTotal  [22m[39m[0m[1mTime[22m
    CHEME 5450 problem set 2 test suite | [32m  12  [39m[36m   12  [39m[0m0.7s
      Setup                             | [32m   7  [39m[36m    7  [39m[0m0.7s
      Calculation                       | [32m   2  [39m[36m    2  [39m[0m0.0s
      Discussion questions              | [32m   3  [39m[36m    3  [39m[0m0.0s

