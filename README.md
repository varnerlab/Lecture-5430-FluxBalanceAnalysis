### Lecture-5430-FluxBalanceAnalysis
This repository holds flux balance analysis lecture for CHEME 5430. 
The lecture is structured in a [Pluto](https://github.com/fonsp/Pluto.jl) notebook and use the [Julia](https://julialang.org) programming language. 

### Installing Julia and Pluto
[Julia](https://julialang.org) is open source, free and runs on all major operating systems and platforms. To install 
[Julia](https://julialang.org) and [Pluto](https://github.com/fonsp/Pluto.jl) please check out the tutorial for 
[MIT 18.S191/6.S083/22.S092 course from Fall 2020](https://computationalthinking.mit.edu/Fall20/installation/).

1. [Install Julia (we are using v1.6.x, newer versions of Julia should also work)](https://julialang.org/downloads/)
1. [Install Pluto.jl](https://github.com/fonsp/Pluto.jl#installation)
1. Clone this repo:
    ```bash
    git clone https://github.com/varnerlab/ENGRI-1120-Cornell-Varner.git
    ```
1. From the Julia REPL (`julia`), run Pluto (a web browser window will pop-up):
    ```julia
    julia> using Pluto
    julia> Pluto.run()
    ```
    _Or you can simply the following in a terminal:_
    ```bash
    julia -E "using Pluto; Pluto.run()"
    ```
2. From Pluto, open the `ConstraintBasedMethods.jl` notebook file located in the root directory (and then enjoy)!

### Alternatively:
* You can view a static `html` version of the lecture [here](https://htmlview.glitch.me/?https://github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis/blob/main/exports/CHEME-5430-ConstraintBasedMethods-Fall-2023.html)
* You can also view a static `pdf` version of the lecture [here](exports/CHEME-5430-ConstraintBasedMethods-Fall-2023.pdf)
