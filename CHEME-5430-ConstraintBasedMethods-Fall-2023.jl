### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ f4ee0417-a667-4594-9c25-a27abebd056e
begin
	
	# Setup the Julia environment -
	using PlutoUI
	using MAT
	
	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_DATA = joinpath(_PATH_TO_NOTEBOOK,"data")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK,"figs")
	
	# return -
	nothing
end

# ╔═╡ 276d790c-99dc-404b-b700-ac36be520aeb
md"""
## Introduction to Metabolic Engineering

Metabolic engineering involves optimizing genetic and regulatory processes within cells or cell-free networks to enhance the production of a desired molecule or protein. This is achieved by manipulating the biochemical networks that convert raw materials into product molecules. The primary aim of metabolic engineering is to:

1. Mathematically model biochemical networks, calculate the yield (product divided substrate) of useful products and identify parts of the network that constrain the production of the products of interest. 
1. Use genetic engineering techniques to modify the biochemical network to relieve constraints limiting production. The modified network can then be modeled to calculate the new product yield and to identify new constraints (back to 1).

Resources for biochemical network information:
* [Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 2000 Jan 1;28(1):27-30. doi: 10.1093/nar/28.1.27. PMID: 10592173; PMCID: PMC102409.](https://www.genome.jp/kegg/)

* [Karp, Peter D et al. “The BioCyc collection of microbial genomes and metabolic pathways.” Briefings in bioinformatics vol. 20,4 (2019): 1085-1093. doi:10.1093/bib/bbx085](https://pubmed.ncbi.nlm.nih.gov/29447345/)

* [Gama-Castro, Socorro, et al. “RegulonDB version 9.0: high-level integration of gene regulation, coexpression, motif clustering and beyond.” Nucleic acids research vol. 44, D1 (2016): D133-43. doi:10.1093/nar/gkv1156](https://pubmed.ncbi.nlm.nih.gov/26527724/)
"""

# ╔═╡ 3f47e0bd-c201-4c82-a392-f912ef36dc3e
md"""
__Fig 1.__ The overall metabolic map from the KEGG database. Each dot (_node_) is a metabolite, each line (_edge_) is a metabolic reaction. 
"""

# ╔═╡ e8b8d393-44b1-4146-b6d5-ec97c6152029
PlutoUI.LocalResource("./figs/KEGG-map01100.png")

# ╔═╡ 0d51af01-c4fb-4f6e-97e1-7aeaf5ba3dec
md"""
## What is Flux Balance Analysis (FBA)?
Flux balance analysis (FBA) is a mathematical modeling and analysis approach which computes the flow (or _flux_) of carbon and energy throughout a metabolic network. FBA, a member of the constraint based family of mathematical modeling tools, is a widely used approach to compute metabolic flux. However, there are alternatives to FBA, such as metabolic flux analysis (MFA), but these alternatives vary more in the solution approach than the structure of the estimation problem. 

Let's look at the following reference to better understand the different components of a flux balance anaysis problem:

* [Orth, J., Thiele, I. & Palsson, B. What is flux balance analysis?. Nat Biotechnol 28, 245–248 (2010). https://doi.org/10.1038/nbt.1614](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3108565/)

Computational tools for working with flux balance analysis models:

* [Heirendt, Laurent et al. “Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0.” Nature protocols vol. 14,3 (2019): 639-702. doi:10.1038/s41596-018-0098-2](https://pubmed.ncbi.nlm.nih.gov/30787451/)

#### Flux balance analysis problem structure
The FBA problem is typically encoded as a linear programming (LP) problem of the form:

$$\max_{v}\sum_{i=1}^{\mathcal{R}}c_{i}v_{i}$$

subject to the constraints:

$$\begin{eqnarray}
\mathbf{S}\mathbf{v} &=& \mathbf{0}\\
\mathcal{L}\leq&\mathbf{v}&\leq\mathcal{U}\\
\dots &=& \dots
\end{eqnarray}$$

where $\mathbf{S}$ denotes the stoichiometric matrix, $c_{i}$ denote the objective coefficients, 
$\mathbf{v}$ denotes the metabolic flux (the unknown that we are trying to estimate), and 
$\mathcal{L}$ ($\mathcal{U}$) denote the permissible lower (upper) bounds on the _unkown_ metabolic flux. The first set of constraints enforces conservation of mass, while the second imparts thermodynamic and kinetic information into the calculation. Finally, there are potentially other types of constraints (both linear and nonlinear) that can be used in this type of problem (we will not cover these here, but these additional constraints may be important in specific applications).

### What is a stoichiometric matrix?
The basis for flux balance calculations is the stoichiometric matrix. 
The stochiometric matrix is a digitial representation of the biochemistry that can occur insides cells (such as that shown in Fig. 1). The stochiometric matrix, denoted by $\mathbf{S}$, is a $\mathcal{M}\times\mathcal{R}$ array of stoichiometric coefficients, denoted by the symbol $\sigma_{ij}$ ($i$ is the row index, $j$ is the col index). Each of the $\mathcal{M}$ rows of $\mathbf{S}$ describes a particular metabolite (_node_ in the metabolic network) while each of the $\mathcal{R}$ columns corresponds to a metabolic reaction (_edge_ in the metabolic network). Thus:

* A stoichiometric coefficient $\sigma_{ij}$ > 0 implies that metabolite $i$ is __produced__ by reaction $j$
* A stoichiometric coefficient $\sigma_{ij}$  = 0 implies that metabolite $i$ is __not connected__ to reaction $j$
* A stoichiometric coefficient $\sigma_{ij}$ < 0 implies that metabolite $i$ is __consumed__ by reaction $j$
"""

# ╔═╡ 43d2f9da-4bed-4c0b-877d-a437a1237c78
md"""

__Example__: Core stoichiometric matrix from _Escherichia coli_: [Orth, Jeffrey D et al. “Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide.” EcoSal Plus vol. 4,1 (2010): 10.1128/ecosalplus.10.2.1. doi:10.1128/ecosalplus.10.2.1](https://pubmed.ncbi.nlm.nih.gov/26443778/)
"""

# ╔═╡ 36ea6b91-6389-45dd-899c-24bfeb265686
begin
	
	# what is the model name -
	model_file_name = "modelReg.mat"
	model_name = "modelReg"
	
	# where is model file?
	_PATH_TO_MODEL_FILE = joinpath(_PATH_TO_DATA, model_file_name)
	
	# load the mat file -> get the cobra_dictionary
	file = matopen(_PATH_TO_MODEL_FILE)
    cobra_dictionary = read(file, model_name)
    close(file)
end

# ╔═╡ 5b3a1883-b7d7-4ba1-82bd-f2db85e2b5fc
keys(cobra_dictionary)

# ╔═╡ 4fda3fa1-166c-4680-9fff-627bb41da9e8
# get the stoichiometric matrix 
stm_sparse = cobra_dictionary["S"]

# ╔═╡ 61d6d920-41f3-49d9-a27e-2ec3f99d2ad8
stm_full = Matrix(stm_sparse)

# ╔═╡ 5f2c354b-af15-495f-ad74-7751bf69c124
mets = cobra_dictionary["metFormulas"]

# ╔═╡ 2b74cc46-2fd3-48fa-b081-e513fae8290a
rxns = cobra_dictionary["rxns"]

# ╔═╡ 003236e7-7f33-4280-b0db-458ba0cebb36
md"""
Additional resources for the stoichiometric matrix

* [Lecture 3: The genome reconstruction process. Systems Biology: Constraint-based Reconstruction and Analysis, Cambridge University Press, 2015](https://www.youtube.com/watch?v=jcRVSJj3XP8&feature=youtu.be)

* [Lecture 9: Properties of the Stoichiometric Matrix. Systems Biology: Constraint-based Reconstruction and Analysis, Cambridge University Press, 2015](https://www.youtube.com/watch?v=0nsjKx6a7X4&feature=youtu.be)

"""

# ╔═╡ 6374e176-5b5b-4dfd-9e12-585fc7317ffb
md"""
### What are flux bounds constraints?
Flux bounds constraints control the permissible ranges that a metabolic reaction rate (flux) can have. Flux bounds constraints can incorporate thermodynamic information (whether a reaction is reversible) as well as kinetic information. 

Resources for thermodynamics in the contraint based world: 

* [Peres, Sabine, and Vincent Fromion. “Thermodynamic Approaches in Flux Analysis.” Methods in molecular biology (Clifton, N.J.) vol. 2088 (2020): 359-367. doi:10.1007/978-1-0716-0159-4_17](https://pubmed.ncbi.nlm.nih.gov/31893383/)

* [Flamholz, Avi et al. “eQuilibrator--the biochemical thermodynamics calculator.” Nucleic acids research vol. 40,Database issue (2012): D770-5. doi:10.1093/nar/gkr874](https://pubmed.ncbi.nlm.nih.gov/22064852/)

Let's do some sample calculations using [eQuilibrator](https://equilibrator.weizmann.ac.il)

Thermodynamic information tells you the directionality of the reaction (is the reaction reversible), but to get the permissible magnitudes of the reaction we need to compute the kinetics. We (and others) formulate flux bounds as the product of succesive corrections t the maximum possible rate e.g., for an irrereversible reaction:

$$0\leq{v_{i}}\leq{V_{max,j}^{\circ}}\left(\frac{\epsilon}{\epsilon^{\circ}}\right)\theta_{j}\left(\dots\right){f_{j}\left(\dots\right)}$$

where $V_{max,j}^{\circ}$ denotes the maximum reaction velocity computed for enzyme $j$ at some characteristic enzyme concetration (and full activity), the ratio $\epsilon/\epsilon^{\circ}$ is a correction for enzyme concentration, $\theta_{j}\left(\dots\right)\in\left[0,1\right]$ is an enzyme activity function (or measurement) and $f_{j}\left(\dots\right)$ is a function describing the substrate dependence of the reaction rate $j$. Both $\theta_{j}\left(\dots\right)$ and $f_{j}\left(\dots\right)$ could have associated parameters, e.g., saturation or binding constants, etc.  

"""

# ╔═╡ 266948f9-b15c-47ae-823c-324e2d9b9573
md"""
### What is the objective of an _E. coli_ cell?

__Short answer__: maximize the growth rate  

__Longer answer__: This question was raging in the metabolic engineering community twenty years ago, but is mostly forgotten today. First, there are long standing arguments from ecology (or even from the perspective of microeconomics) as to why _E.coli_ would be goal oriented. However, Palsson and coworkers largely put this question to rest (at least for _E.coli_) with the publication:

* [Ibarra, Rafael U et al. “Escherichia coli K-12 undergoes adaptive evolution to achieve in silico predicted optimal growth.” Nature vol. 420,6912 (2002): 186-9. doi:10.1038/nature01149](https://pubmed.ncbi.nlm.nih.gov/12432395/)

__What about cell-free networks__? 
No cells, but most of the biochemistry. Is there an objective? This is an open question. However, initial work in this area suggests the maximization of translation, at least for single gene systems:

* [Vilkhovoy, Michael et al. “Sequence Specific Modeling of E. coli Cell-Free Protein Synthesis.” ACS synthetic biology vol. 7,8 (2018): 1844-1857. doi:10.1021/acssynbio.7b00465](https://pubmed.ncbi.nlm.nih.gov/29944340/)

"""

# ╔═╡ 3518d9b8-af6f-4a6b-8a11-74d5bf9dbbde
md"""
## Flux balance analysis case studies and examples

* [Wayman, Joseph A et al. “Improving designer glycan production in Escherichia coli through model-guided metabolic engineering.” Metabolic engineering communications vol. 9 e00088. 29 Mar. 2019, doi:10.1016/j.mec.2019.e00088](https://pubmed.ncbi.nlm.nih.gov/31008057/)

"""

# ╔═╡ 01bb1ad0-616b-11ee-0cdf-af81f2517005
html"""
<style>
main {
    max-width: 1200px;
  	width: 85%;
	margin: auto;
	font-family: "Roboto, monospace";
}

a {
	color: blue;
	text-decoration: none;
}

.H1 {
    padding: 0px 30px;
}
</style>
<script>
var section = 0;
var subsection = 0;
var headers = document.querySelectorAll('h2, h3');
for (var i=0; i < headers.length; i++) {
    var header = headers[i];
    var text = header.innerText;
    var original = header.getAttribute("text-original");
    if (original === null) {
        // Save original header text
        header.setAttribute("text-original", text);
    } else {
        // Replace with original text before adding section number
        text = header.getAttribute("text-original");
    }
    var numbering = "";
    switch (header.tagName) {
        case 'H2':
            section += 1;
            numbering = section + ".";
            subsection = 0;
            break;
        case 'H3':
            subsection += 1;
            numbering = section + "." + subsection;
            break;
    }
    header.innerText = numbering + " " + text;
};
</script>"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
MAT = "~0.10.6"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "1c1565d1bd164b6de99b3edd13c653fc591c3e98"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BufferedStreams]]
git-tree-sha1 = "4ae47f9a4b1dc19897d3743ff13685925c5202ec"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.2.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "26407bd1c60129062cec9da63dc7d08251544d53"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.1"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "38c8874692d48d5440d5752d6c74b0c6b0b60739"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.2+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "ed1cf0a322d78cee07718bed5fd945e2218c35a1"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.6"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "8a5b4d2220377d1ece13f49438d71ad20cf1ba83"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.1.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "781916a2ebf2841467cda03b6f1af43e23839d85"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.9"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "6979eccb6a9edbbb62681e158443e79ecc0d056a"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.3.1+0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a7023883872e52bc29bcaac74f19adf39347d2d5"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "f3080f4212a8ba2ceb10a34b938601b862094314"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.5+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ceeda72c9fd6bbebc4f4f598560789145a8b6c4c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.11+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eddd19a8dea6b139ea97bdc8a0e2667d4b661720"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.0.6+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─276d790c-99dc-404b-b700-ac36be520aeb
# ╟─3f47e0bd-c201-4c82-a392-f912ef36dc3e
# ╟─e8b8d393-44b1-4146-b6d5-ec97c6152029
# ╟─0d51af01-c4fb-4f6e-97e1-7aeaf5ba3dec
# ╟─43d2f9da-4bed-4c0b-877d-a437a1237c78
# ╠═36ea6b91-6389-45dd-899c-24bfeb265686
# ╠═5b3a1883-b7d7-4ba1-82bd-f2db85e2b5fc
# ╠═4fda3fa1-166c-4680-9fff-627bb41da9e8
# ╠═61d6d920-41f3-49d9-a27e-2ec3f99d2ad8
# ╠═5f2c354b-af15-495f-ad74-7751bf69c124
# ╠═2b74cc46-2fd3-48fa-b081-e513fae8290a
# ╟─003236e7-7f33-4280-b0db-458ba0cebb36
# ╟─6374e176-5b5b-4dfd-9e12-585fc7317ffb
# ╟─266948f9-b15c-47ae-823c-324e2d9b9573
# ╟─3518d9b8-af6f-4a6b-8a11-74d5bf9dbbde
# ╟─f4ee0417-a667-4594-9c25-a27abebd056e
# ╟─01bb1ad0-616b-11ee-0cdf-af81f2517005
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
