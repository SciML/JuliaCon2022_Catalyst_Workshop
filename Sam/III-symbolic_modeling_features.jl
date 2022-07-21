### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 3177b2e3-eefd-4faf-b292-2c7c6f0bbe22
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, ModelingToolkit, DifferentialEquations, Plots
end

# ╔═╡ bcdd2649-93c4-4d20-a4f6-e4437b8b842b
using Distributions: Geometric

# ╔═╡ 192af62c-08fc-11ed-35a9-bb2a5d18c32f
html"<button onclick=present()>Present</button>"

# ╔═╡ 7cb53ff0-1751-4940-a4b3-bec68e81a005
html"""
<script>
    const calculate_slide_positions = (/** @type {Event} */ e) => {
        const notebook_node = /** @type {HTMLElement?} */ (e.target)?.closest("pluto-editor")?.querySelector("pluto-notebook")
		console.log(e.target)
        if (!notebook_node) return []
        const height = window.innerHeight
        const headers = Array.from(notebook_node.querySelectorAll("pluto-output h1, pluto-output h2"))
        const pos = headers.map((el) => el.getBoundingClientRect())
        const edges = pos.map((rect) => rect.top + window.pageYOffset)
        edges.push(notebook_node.getBoundingClientRect().bottom + window.pageYOffset)
        const scrollPositions = headers.map((el, i) => {
            if (el.tagName == "H1") {
                // center vertically
                const slideHeight = edges[i + 1] - edges[i] - height
                return edges[i] - Math.max(0, (height - slideHeight) / 2)
            } else {
                // align to top
                return edges[i] - 20
            }
        })
        return scrollPositions
    }
    const go_previous_slide = (/** @type {Event} */ e) => {
        const positions = calculate_slide_positions(e)
        const pos = positions.reverse().find((y) => y < window.pageYOffset - 10)
        if (pos) window.scrollTo(window.pageXOffset, pos)
    }
    const go_next_slide = (/** @type {Event} */ e) => {
        const positions = calculate_slide_positions(e)
        const pos = positions.find((y) => y - 10 > window.pageYOffset)
        if (pos) window.scrollTo(window.pageXOffset, pos)
    }
	const left_button = document.querySelector(".changeslide.prev")
	const right_button = document.querySelector(".changeslide.next")
	left_button.addEventListener("click", go_previous_slide)
	right_button.addEventListener("click", go_next_slide)
</script>
"""

# ╔═╡ 32897144-3033-45bb-a878-f72f6882df20
md"
```math
\DeclareMathOperator*{\argmin}{arg\,min}
\DeclareMathOperator*{\rank}{rank}
\DeclareMathOperator*{\col}{col}
\DeclareMathOperator*{\proj}{proj}
\DeclareMathOperator*{\sign}{sign}
\DeclareMathOperator*{\span}{span}
\DeclareMathOperator*{\cond}{cond}
\DeclareMathOperator*{\repart}{Re}
\DeclareMathOperator*{\impart}{Im}
\newcommand{\vec}[1]{\boldsymbol{#1}}
\newcommand{\vx}{\vec{x}}
\newcommand{\vy}{\vec{y}}
\newcommand{\vz}{\vec{z}}
\newcommand{\vb}{\vec{b}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\D}[2]{\frac{d#1}{d#2}}
\newcommand{\PD}[2]{\frac{\partial#1}{\partial#2}}
\newcommand{\PDD}[2]{\frac{\partial^{2}{#1}}{\partial{#2}^{2}}}
\newcommand{\paren}[1]{\left(#1\right)}
```
"

# ╔═╡ 146fee58-6588-4bfc-94d9-cfa7302b59a0
md"Package setup:"

# ╔═╡ dbfd836f-d475-48e7-9afd-8f02d9785f42
md"## _*Symbolic Modeling Features of Catalyst.jl*_
Catalyst offers a lot of additional functionality that can be exploited by using the underlying ModelingToolkit.jl and Symbolics.jl CAS features:
- Building networks symbolically and programmatically
    - Interpolation into the DSL
    - `@reaction` macro
- Parametric stoichiometry
- External Julia functions
- Compositional modeling
- Adding constraint equations
"

# ╔═╡ e63ce294-4e50-4fec-8092-8f34facfd7f5
md"## _*Symbolic Interface*_
Let's see how we can define the Oregonator model via the symbolic interface. First we specify it via the DSL approach we previously learned:"

# ╔═╡ 6c8e9c23-eca9-4b91-8391-8eefeb3b3745
rn = @reaction_network oregonator begin
	k₁, X + Y --> ∅
	k₂, Y --> X
	k₃, X --> 2X + Z
	k₄, 2X --> ∅
	k₅, Z --> Y
end k₁ k₂ k₃ k₄ k₅

# ╔═╡ 0c4f8cd6-b05e-4386-956f-9eaa404ab656
md"Let's first take a quick look at what the ODE solutions look like:"

# ╔═╡ 1cb68f8a-e33e-4c67-90d6-94f51e6dafab
let
	u₀ = [:X => 500, :Y => 1000, :Z => 2100]
	p = [:k₁ => .1, :k₂ => 2.0, :k₃ => 104.0, :k₄ => .016, :k₅ => 26.0]
	tspan = (0.0, 5.0)
	oprob = ODEProblem(rn, u₀, tspan, p)
	sol = solve(oprob, Tsit5())
	plot(sol)
end

# ╔═╡ 7cd8a85f-d5fd-45dd-8254-941d7b6afabd
md"## _*Building the Oregonator Symbolically*_
We first define our symbolic parameters:"

# ╔═╡ 47ca3da2-28db-47dd-98e8-7f5bec0c9597
@parameters k₁ k₂ k₃ k₄ k₅

# ╔═╡ bbcf38ea-518c-4365-82a1-7d9eed1536a4
md"We next define our state variables, the species populations"

# ╔═╡ b042d589-0c59-4243-b6d9-1e4a55dbcc43
@variables t X(t) Y(t) Z(t)

# ╔═╡ d7184781-18d3-4c53-b9a6-8ab628a4d603
md"We have now constructed [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) symbolic parameters and variables we can manipulate in Julia using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl):"

# ╔═╡ 8ff33d17-902c-44bb-ba24-9e89c06fb887
ex = X^2 + sin(Y) / Z

# ╔═╡ d9393282-8389-43af-8ab9-f9d9f3943cad
simplify((X-1)^2 / (X-1))

# ╔═╡ a052955b-e291-4618-9ff5-df9dd630d753
md"Now we can construct each reaction in the Oregonator symbolically in two ways. Let's recall our reactions:"

# ╔═╡ 53c6cf2d-3b56-460b-93ec-1b930e1e89d3
rn

# ╔═╡ 82d862db-0cf0-4cc5-973d-2e832dcb39c7
md"Symbolic method 1:"

# ╔═╡ e1266df7-84bf-48d5-ac2d-e1dec70ed7e4
rxs = [Reaction(k₁, [X, Y], nothing),
	   Reaction(k₂, [Y], [X]),
	   Reaction(k₃, [X], [X, Z], [1], [2,1]),
	   Reaction(k₄, [X], nothing, [2], nothing),
	   Reaction(k₅, [Z], [Y])]

# ╔═╡ 1321fde4-c779-4720-9fc1-803503f18bab
md"Notice the general format is
```julia
Reaction(rate_constant_expression, substrates, products, substoichs, prodstoichs)
```
with `nothing` used where there are no substrates and/or products.

Catalyst also provides a simple macro that can be used to write reactions in the DSL notation one at a time:
"

# ╔═╡ 2799b012-76c2-474e-ac57-d11f0320c7fc
rxs2 = [(@reaction k₁, X + Y --> ∅),
	    (@reaction k₂, Y --> X),
	    (@reaction k₃, X --> 2X + Z),
	    (@reaction k₄, 2X --> ∅),
	    (@reaction k₅, Z --> Y)]

# ╔═╡ 250980aa-9e10-4b74-93f4-f408c3740d8d
md"Let's check the two reaction lists are equivalent:"

# ╔═╡ a1cb0d70-784a-4253-acf0-de1bbd720ec6
rxs == rxs2

# ╔═╡ a8e596a6-ef9f-45b1-b47b-3b99eef5bdac
md"Now we are ready to construct our `ReactionSystem`, which stores the chemical reaction network and is the output of the `@reaction_network` macro:"

# ╔═╡ 68be114a-471e-44a8-a0b0-7ebc8b2c96c8
@named oregonator = ReactionSystem(rxs, t)

# ╔═╡ f21d53e9-1674-47db-9b63-d78c9dfedf2e
oregonator == rn

# ╔═╡ 2cbc666d-1b4a-49f7-9e0e-b1a8fc9ca18f
md"*Note, equality between `ReactionSystem`s requires they have the same name, and looking above we see we did give `rn` the name `oregonator` too:*"

# ╔═╡ 57f7d7ce-0b61-424a-a4cd-a6750a45cc8d
nameof(rn), nameof(oregonator)

# ╔═╡ b130ca6b-deb1-475c-983b-5e1c6bcf6ab1
md"## _*Parametric Stoichiometry*_
So far we have simply used whole numbers to represent stoichiometry in reactions, however, Catalyst actually allows non-negative real numbers and even more general expressions. 

Let's look at a model with \"bursty\" production of proteins from a gene
1. We have a gene that can flip between an off state, ``G_-``, and an on state, ``G_+``.
2. Binding of protein dimers is what flips the gene between the two states.
3. Proteins can degrate away.
4. In the on state, the gene expression can occur, and leads to the production of a random amount of protein, distributed accoding to a (shifted) geometric random variable. 

Let `b` be our parameter that sets the average burst size. Our parameters and variables are then:
"

# ╔═╡ 60dfb1e3-ef92-465e-b1d7-a351036323a9
@parameters b k₊ k₋ kₚ γₚ

# ╔═╡ 73c18851-f5e1-4a17-ae19-0fef083ed5e8
@variables G₋(t) G₊(t) P(t)

# ╔═╡ b84c5bff-c4e3-41ed-8686-6253b4b87673
md"Let's now define our stochiometric coefficient, `m`, that represents the random amount of protein produced during a transcription event. We use [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) to sample this variable:"

# ╔═╡ 1f353a66-436d-4ffc-933d-ca3208c7c391
md"First, to make Distributions' geometric distribution accessible in Symbolics.jl expression we need to register the function, which tells Symbolics not to try to evaluate this function with symbolic variables"

# ╔═╡ 56c15ca4-1a26-4f7f-9d4a-dbe138614f03
@register_symbolic Geometric(b)

# ╔═╡ b4ad6737-55d8-417c-bc31-04898807089b
m = rand(Geometric(1/b)) + 1

# ╔═╡ 7f6e982b-8710-4d1a-84ee-0ed77e568a8c
md"This corresponds to a symbolic variable that when evaluated generates a sample from the Geometric Distribution with mean `b`, note we shift by one to obtain a sample from the Shifted Geometric Distribution!"

# ╔═╡ 306bcb0c-5134-4e50-a72a-5d559a1726d7


# ╔═╡ 63fd312f-50da-46e5-934a-6deacfbd474a


# ╔═╡ e7a89226-4757-42c1-825a-b117ae439c1e


# ╔═╡ 7a065cca-ffac-47cb-bc29-c00237d6936a
md"## _*Appendix*_"

# ╔═╡ 0bbeb157-ebc7-4c99-aa1d-0adbeeb0f27f
# plot defaults
begin
	_fsz = 12
	_tsz = 12
	default(size=(800,400), 
			xtickfontsize=_tsz,
			ytickfontsize=_tsz,
		    titlefontsize=_fsz,
	     	xguidefontsize=_fsz, 
		    yguidefontsize=_fsz,
		    legendfontsize=_fsz,
	        margin=5Plots.mm);
end

# ╔═╡ Cell order:
# ╟─192af62c-08fc-11ed-35a9-bb2a5d18c32f
# ╟─7cb53ff0-1751-4940-a4b3-bec68e81a005
# ╟─32897144-3033-45bb-a878-f72f6882df20
# ╟─146fee58-6588-4bfc-94d9-cfa7302b59a0
# ╠═3177b2e3-eefd-4faf-b292-2c7c6f0bbe22
# ╟─dbfd836f-d475-48e7-9afd-8f02d9785f42
# ╟─e63ce294-4e50-4fec-8092-8f34facfd7f5
# ╠═6c8e9c23-eca9-4b91-8391-8eefeb3b3745
# ╟─0c4f8cd6-b05e-4386-956f-9eaa404ab656
# ╠═1cb68f8a-e33e-4c67-90d6-94f51e6dafab
# ╟─7cd8a85f-d5fd-45dd-8254-941d7b6afabd
# ╠═47ca3da2-28db-47dd-98e8-7f5bec0c9597
# ╟─bbcf38ea-518c-4365-82a1-7d9eed1536a4
# ╠═b042d589-0c59-4243-b6d9-1e4a55dbcc43
# ╟─d7184781-18d3-4c53-b9a6-8ab628a4d603
# ╠═8ff33d17-902c-44bb-ba24-9e89c06fb887
# ╠═d9393282-8389-43af-8ab9-f9d9f3943cad
# ╟─a052955b-e291-4618-9ff5-df9dd630d753
# ╟─53c6cf2d-3b56-460b-93ec-1b930e1e89d3
# ╟─82d862db-0cf0-4cc5-973d-2e832dcb39c7
# ╠═e1266df7-84bf-48d5-ac2d-e1dec70ed7e4
# ╟─1321fde4-c779-4720-9fc1-803503f18bab
# ╠═2799b012-76c2-474e-ac57-d11f0320c7fc
# ╟─250980aa-9e10-4b74-93f4-f408c3740d8d
# ╠═a1cb0d70-784a-4253-acf0-de1bbd720ec6
# ╟─a8e596a6-ef9f-45b1-b47b-3b99eef5bdac
# ╠═68be114a-471e-44a8-a0b0-7ebc8b2c96c8
# ╠═f21d53e9-1674-47db-9b63-d78c9dfedf2e
# ╟─2cbc666d-1b4a-49f7-9e0e-b1a8fc9ca18f
# ╠═57f7d7ce-0b61-424a-a4cd-a6750a45cc8d
# ╟─b130ca6b-deb1-475c-983b-5e1c6bcf6ab1
# ╠═60dfb1e3-ef92-465e-b1d7-a351036323a9
# ╠═73c18851-f5e1-4a17-ae19-0fef083ed5e8
# ╠═b84c5bff-c4e3-41ed-8686-6253b4b87673
# ╠═bcdd2649-93c4-4d20-a4f6-e4437b8b842b
# ╠═1f353a66-436d-4ffc-933d-ca3208c7c391
# ╠═56c15ca4-1a26-4f7f-9d4a-dbe138614f03
# ╠═b4ad6737-55d8-417c-bc31-04898807089b
# ╠═7f6e982b-8710-4d1a-84ee-0ed77e568a8c
# ╠═306bcb0c-5134-4e50-a72a-5d559a1726d7
# ╠═63fd312f-50da-46e5-934a-6deacfbd474a
# ╠═e7a89226-4757-42c1-825a-b117ae439c1e
# ╟─7a065cca-ffac-47cb-bc29-c00237d6936a
# ╠═0bbeb157-ebc7-4c99-aa1d-0adbeeb0f27f
