### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 45376f68-2552-11ec-32d7-a5de0d6b6882
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

# ╔═╡ 16a831a8-14b2-4d1e-86e7-a28f21d0dbd4
md"""
## Introduction to Metabolic Engineering

Metabolic engineering is the practice of optimizing genetic and regulatory processes within cells to increase the cell's production of a desired molecule or protein of interest. 

Metabolic engineers manipulate the biochemical networks used by cells to convert raw materials into molecules necessary for the cell's survival. Metabolic engineering specifically seeks to:

1. Mathematically model biochemical networks, calculate the yield (product divided substrate) of useful products, and identify parts of the network that constrain the production of the products of interest. 
1. Use genetic engineering techniques to modify the biochemical network in order to relieve constraints limiting production. The modified network can then be modeled to calculate the new product yield, and to identify new constraints (back to 1).

Resources for biochemical network information:
* [Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 2000 Jan 1;28(1):27-30. doi: 10.1093/nar/28.1.27. PMID: 10592173; PMCID: PMC102409.](https://www.genome.jp/kegg/)

* [Karp, Peter D et al. “The BioCyc collection of microbial genomes and metabolic pathways.” Briefings in bioinformatics vol. 20,4 (2019): 1085-1093. doi:10.1093/bib/bbx085](https://pubmed.ncbi.nlm.nih.gov/29447345/)

* [Gama-Castro, Socorro et al. “RegulonDB version 9.0: high-level integration of gene regulation, coexpression, motif clustering and beyond.” Nucleic acids research vol. 44,D1 (2016): D133-43. doi:10.1093/nar/gkv1156](https://pubmed.ncbi.nlm.nih.gov/26527724/)

"""

# ╔═╡ 716d10ed-38bd-4654-9a54-33dfba020591
md"""
__Fig 1.__ The overall metabolic map from the KEGG database. Each dot (_node_) is a metabolite, each line (_edge_) is a metabolic reaction. 
"""

# ╔═╡ ad2d1e2b-4222-43d0-819a-4a7d14746327
PlutoUI.LocalResource("./figs/KEGG-map01100.png")

# ╔═╡ cef733a0-5541-4338-bd28-0916b697fb19
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

# ╔═╡ 679f15d4-e94d-43eb-8342-701111c46201
md"""

__Example__: Core stoichiometric matrix from _Escherichia coli_: [Orth, Jeffrey D et al. “Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide.” EcoSal Plus vol. 4,1 (2010): 10.1128/ecosalplus.10.2.1. doi:10.1128/ecosalplus.10.2.1](https://pubmed.ncbi.nlm.nih.gov/26443778/)
"""

# ╔═╡ 0f683f10-71d7-477c-9b7c-49545e60fa8a
# load the stoichiometric matrix, and mess around with it -
begin
	
	# what is the model name -
	model_file_name = "modelReg.mat"
	model_name = "modelReg"
	
	# where is model file?
	_PATH_TO_MODEL_FILE = joinpath(_PATH_TO_DATA,model_file_name)
	
	# load the mat file -> get the cobra_dictionary
	file = matopen(_PATH_TO_MODEL_FILE)
    cobra_dictionary = read(file, model_name)
    close(file)
end

# ╔═╡ 095a2f1a-6ae3-451a-bc3f-7dad9a7834ea
keys(cobra_dictionary)

# ╔═╡ 5d7c052a-73d7-487a-939e-f813561daff9
# get the stoichiometric matrix 
stm_sparse = cobra_dictionary["S"]

# ╔═╡ b3e23596-0b67-4ed1-8caa-6353307ee315
stm_full = Matrix(stm_sparse)

# ╔═╡ 5f19c57c-e4ac-4c8c-9bc0-f6a015c8fb3f
mets = cobra_dictionary["metFormulas"]

# ╔═╡ 1062535d-3453-4a00-89f5-66060c150dbf
rxns = cobra_dictionary["rxns"]

# ╔═╡ 5578b9b1-8fbd-4c22-a89d-cd39ea63fe87
md"""
Additional resources for the stoichiometric matrix

* [Lecture 3: The genome reconstruction process. Systems Biology: Constraint-based Reconstruction and Analysis, Cambridge University Press, 2015](https://www.youtube.com/watch?v=jcRVSJj3XP8&feature=youtu.be)

* [Lecture 9: Properties of the Stoichiometric Matrix. Systems Biology: Constraint-based Reconstruction and Analysis, Cambridge University Press, 2015](https://www.youtube.com/watch?v=0nsjKx6a7X4&feature=youtu.be)

"""

# ╔═╡ 6629d6a5-af02-40c4-a7a3-a7146bf6638b
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

# ╔═╡ baedaf07-7236-432c-a67b-bd7eaa99bb7c
md"""
### What is the objective of an _E. coli_ cell?

__Short answer__: maximize the growth rate  

__Longer answer__: This question was raging in the metabolic engineering community twenty years ago, but is mostly forgotten today. First, there are long standing arguments from ecology (or even from the perspective of microeconomics) as to why _E.coli_ would be goal oriented. However, Palsson and coworkers largely put this question to rest (at least for _E.coli_) with the publication:

* [Ibarra, Rafael U et al. “Escherichia coli K-12 undergoes adaptive evolution to achieve in silico predicted optimal growth.” Nature vol. 420,6912 (2002): 186-9. doi:10.1038/nature01149](https://pubmed.ncbi.nlm.nih.gov/12432395/)

__What about cell-free networks__? 
No cells, but most of the biochemistry. Is there an objective? This is an open question. However, initial work in this area suggests the maximization of translation, at least for single gene systems:

* [Vilkhovoy, Michael et al. “Sequence Specific Modeling of E. coli Cell-Free Protein Synthesis.” ACS synthetic biology vol. 7,8 (2018): 1844-1857. doi:10.1021/acssynbio.7b00465](https://pubmed.ncbi.nlm.nih.gov/29944340/)

"""

# ╔═╡ fd376fb1-1617-48a8-877b-31aa8f85476a
md"""
## Flux balance analysis case studies and examples

* [Wayman, Joseph A et al. “Improving designer glycan production in Escherichia coli through model-guided metabolic engineering.” Metabolic engineering communications vol. 9 e00088. 29 Mar. 2019, doi:10.1016/j.mec.2019.e00088](https://pubmed.ncbi.nlm.nih.gov/31008057/)

"""

# ╔═╡ 774dd760-cebf-44ff-b753-e54353c1317e
TableOfContents(title="CHEME 5430: Constraint Based Methods and Models", depth=4)

# ╔═╡ debd1859-15ce-49ec-b45a-4707695f5c8a
html"""
<style>
main {
    max-width: 1200px;
    margin: auto;
	width: 90%;
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
"""

# ╔═╡ 8904e3f1-88ad-424e-8ac2-e4c8dfd01c4a
html"""
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
</script>
"""

# ╔═╡ d1bd21e8-86c0-4607-953b-d33a7b964ecc
html"""
<style>
main {
    max-width: 1200px;
  	width: 45%;
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
</script>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
MAT = "~0.10.1"
PlutoUI = "~0.7.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "84cf7d0f8fd46ca6f1b3e0305b4b4a37afe50fd6"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.0"

[[Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "e747dac84f39c62aff6956651ec359686490134e"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.0+0"

[[BufferedStreams]]
deps = ["Compat", "Test"]
git-tree-sha1 = "5d55b9486590fdda5905c275bb21ce1f0754020f"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.0.0"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[HDF5]]
deps = ["Blosc", "Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "83173193dc242ce4b037f0263a7cc45afb5a0b85"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.15.6"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "5c62992f3d46b8dce69bdd234279bb5a369db7d5"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.1"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "9d8c00ef7a8d110787ff6f170579846f776133a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.4"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─16a831a8-14b2-4d1e-86e7-a28f21d0dbd4
# ╟─716d10ed-38bd-4654-9a54-33dfba020591
# ╠═ad2d1e2b-4222-43d0-819a-4a7d14746327
# ╟─cef733a0-5541-4338-bd28-0916b697fb19
# ╟─679f15d4-e94d-43eb-8342-701111c46201
# ╠═0f683f10-71d7-477c-9b7c-49545e60fa8a
# ╠═095a2f1a-6ae3-451a-bc3f-7dad9a7834ea
# ╠═5d7c052a-73d7-487a-939e-f813561daff9
# ╠═b3e23596-0b67-4ed1-8caa-6353307ee315
# ╠═5f19c57c-e4ac-4c8c-9bc0-f6a015c8fb3f
# ╠═1062535d-3453-4a00-89f5-66060c150dbf
# ╟─5578b9b1-8fbd-4c22-a89d-cd39ea63fe87
# ╟─6629d6a5-af02-40c4-a7a3-a7146bf6638b
# ╟─baedaf07-7236-432c-a67b-bd7eaa99bb7c
# ╟─fd376fb1-1617-48a8-877b-31aa8f85476a
# ╠═45376f68-2552-11ec-32d7-a5de0d6b6882
# ╠═774dd760-cebf-44ff-b753-e54353c1317e
# ╠═debd1859-15ce-49ec-b45a-4707695f5c8a
# ╠═8904e3f1-88ad-424e-8ac2-e4c8dfd01c4a
# ╟─d1bd21e8-86c0-4607-953b-d33a7b964ecc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
