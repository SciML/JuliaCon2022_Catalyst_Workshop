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

	using Catalyst, ModelingToolkit, DifferentialEquations, Plots, GraphRecipes, IfElse, NonlinearSolve
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

# ╔═╡ d35195ef-de7a-4f4c-8a38-ca9b4ec1b742
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
	        margin=5Plots.mm,
	        lw=2);
end

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
ex = (X^2 - 2*X + 1) / (X-1)

# ╔═╡ d9393282-8389-43af-8ab9-f9d9f3943cad
ex2 = simplify(ex)

# ╔═╡ bf439579-ab39-4e3e-9c1d-0d71ac7d8b9e
substitute(ex2, Dict((X-1) => Y))

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

Let's look at a model with \"bursty\" production of proteins from a gene with the reactions
```math
\begin{aligned}
\varnothing &\overset{k}{\to} m P, \\
P &\overset{\gamma}{\to} \varnothing
\end{aligned}
```
- ``m`` is a shifted geometric random variable for the number of proteins produced. 
- ``b`` is our parameter that sets the average burst size when one transcription event occurs. 
Our parameters and variables are then:
"

# ╔═╡ 60dfb1e3-ef92-465e-b1d7-a351036323a9
@parameters b k γ

# ╔═╡ 73c18851-f5e1-4a17-ae19-0fef083ed5e8
@variables P(t)

# ╔═╡ b84c5bff-c4e3-41ed-8686-6253b4b87673
md"Let's now define our stochiometric coefficient, `m`, that represents the random amount of protein produced during a transcription event. We use [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) to sample this variable:"

# ╔═╡ c3be82d6-54f5-4876-8583-4f834417b0bf
md"Let's try to create a symbolic variable that when evaluated will correspond to sampling from this distribution. 

If `b` were a normal Julia `Float64` we might try the following
```julia
m = rand(Geometric(1/b)) + 1
```
Unfortunately this will give the following error:
```julia
ERROR: TypeError: non-boolean (Num) used in boolean context
```

To make Distributions' geometric distribution accessible in Symbolics.jl expressions we need to register the function, which tells Symbolics not to try to evaluate this function with symbolic variables
"

# ╔═╡ 56c15ca4-1a26-4f7f-9d4a-dbe138614f03
@register_symbolic Geometric(b)

# ╔═╡ b4ad6737-55d8-417c-bc31-04898807089b
m = rand(Geometric(1/b)) + 1

# ╔═╡ 7f6e982b-8710-4d1a-84ee-0ed77e568a8c
md"This corresponds to a symbolic variable that when evaluated generates a sample from the Geometric Distribution with mean `b`, note we shift by one to obtain a sample from the Shifted Geometric Distribution!"

# ╔═╡ 306bcb0c-5134-4e50-a72a-5d559a1726d7
burstyrxs = [Reaction(k, nothing, [P], nothing, [m]),
			 Reaction(γ, [P], nothing)]

# ╔═╡ a5277543-582d-4958-80f4-1f43857601a0
@named burstygene = ReactionSystem(burstyrxs, t)

# ╔═╡ 63fd312f-50da-46e5-934a-6deacfbd474a
burstysol = let 
	p = [:k => 10.0, :γ => 1.0, :b => 10]
	u₀ =[:P => 0]
	tspan = (0., 2) 
	dprob = DiscreteProblem(burstygene, u₀, tspan, p)
	jprob = JumpProblem(burstygene, dprob, Direct())
	sol = solve(jprob, SSAStepper())
end

# ╔═╡ cbb87922-da86-494e-8e46-7ecbad29f432
md"Suppose we just want to plot `P`, we can use symbolic indexing with `sol` to get at it. There are several ways we can do this. 

1. We just use `burstygene.P` to access the symbolic variable:"

# ╔═╡ ff6cd1cf-3753-4f93-94e1-dbbaccf97aa8
plot(burstysol, vars=burstygene.P, legend=:bottomright)

# ╔═╡ 1da19c48-9532-4e48-9bd5-ce968e3de493
md"
2. We `@unpack` the variable from the system `burstygene`"

# ╔═╡ edd33cf5-b00e-4b13-bdfe-175cb01f6822
let
	@unpack P = burstygene
	plot(burstysol, vars=P)
end

# ╔═╡ 92c7f19d-9299-488e-87e9-76c68fe26405
md"3. We just create a symbolic variable named `P(t)` from scratch and use it"

# ╔═╡ e7a89226-4757-42c1-825a-b117ae439c1e
let
	@variables P(t)
	plot(burstysol, vars=P)
end

# ╔═╡ 9252661d-6fb8-47fb-805d-895f228df4e2
md"We can use each of these to access the solution at specific times too, for example:"

# ╔═╡ 03ce0b66-6860-4f8b-9f22-ec9ffaec60d1
u = let
	tv = collect(range(0, 1, length=10))
	burstysol(tv, idxs=P)
end

# ╔═╡ 4430dc04-3a31-408c-85d9-1efe36c088d1
u.u   # P at the times in tv

# ╔═╡ 0fe3d571-acf6-462c-a7f3-07e65ac0849f
u.t   # tv

# ╔═╡ b34914dc-bdbb-4ba0-b992-b4295cf0c642
md"## _*Interpolation into the DSL*_
- It can be useful to mix the simple DSL notation with creating some portions of a system symbolically! 
- We can interpolate symbolic variables and expressions into the DSL using the `$` symbol. 
- Let's redo our previous example using interpolation. 
"

# ╔═╡ 60263535-9291-483f-b28e-fc8c26f4019f
burstygene2 = let
	@parameters b
	m = 1 + rand(Geometric(1/b))
	@reaction_network burstygene begin
		k, ∅ --> $m*P
		γ, P --> ∅
	end k γ	
end

# ╔═╡ 121b11c2-aadc-4e73-871a-d6c79c8fa3d9
md"This is the same system as before:"

# ╔═╡ 95e3bf53-d644-480a-9c75-602501657f7e
burstygene == burstygene2

# ╔═╡ 15dede87-f769-4d16-ac76-95bf4d8777cd
md"One can interpolate previously defined symbolic variables, parameters, or Julia variables storing expressions into the DSL."

# ╔═╡ a4bd43d3-2b45-4a64-b17c-ca188d1d7d34
md"## _*Note about registering functions*_
- `@register_symbolic` allows us to use general Julia functions in symbolic expressions.
- Symbolics.jl will trace into the functions we use in symbolic expressions (i.e. try to evaluate them with symbolic variables). 
  - This allows Symbolics to fully caclulate derivatives in ODE Jacobian functions and such, but doesn't work for all Julia functions. 
For example,
"

# ╔═╡ a891e0ce-0cf5-4e5e-9ffb-b8cf257dd750
let
	@parameters b
	f(x) = 1 + x^3
	f(b)
end

# ╔═╡ 82ac4773-4320-4683-85d1-8bdc0f5a77a3
md"works, but"

# ╔═╡ ed06ac37-98bb-4a28-81c7-29885405fa3b
md"fails. If we register it though it will be usable:"

# ╔═╡ 771e0080-a6d4-4de8-a7ec-3bb3949fc57f
let 
	# this is needed to define a version of it for symbolic variables
	import Base: round   
	@parameters b
	@register_symbolic round(x)
	round(b)
end

# ╔═╡ 5535c9cc-bb87-4046-8621-93214829a012
let
	@parameters b
	round(b)
end

# ╔═╡ 09d51be7-0972-4701-bfd7-471b207f518d
md"
- When registering a function, Symbolics will no longer \"trace\" into it. 
- However, you need to also register a derivative for that function if you want to differentiate expressions that involve it with respect to its arguments!
- This doesn't matter for jump process models, but can arise when calculating Jacobians of ODE models or in various types of analysis (like sensitivity analysis).
"

# ╔═╡ b6dcccae-1aa3-4778-9b5b-bd0d8858df95
md"## _*Compositional Modeling*_
Catalyst supports several methods for composing models together based on the ModelingToolkit API, these include:
- `extend` - which merges several networks together into one network
- `compose` - which creates a higherarchical tree of networks using subsystems.

Let's see how each works
"

# ╔═╡ fac7b265-4802-4fa2-87b8-a08a57e6ec6c
basern = @reaction_network rn1 begin
           k, A + B --> C
         end k

# ╔═╡ becb9021-3e24-47a3-8893-3e5a55242958
newrn = @reaction_network rn2 begin
        r, C --> A + B
      end r

# ╔═╡ 489b6939-5505-4bf7-8e5f-1a48ae4e6fe4
@named extended_network = extend(newrn, basern)

# ╔═╡ 09f90ebf-4fa8-4cca-9fb5-dcfc64084d69
md"
- Here we have merged the two networks into a new network with the same sets of reactions.
- `compose` works hierarchically instead, treating the system's species and parameters as distinct. 

Let's make one more network to help us see how this works
"

# ╔═╡ 0d736c32-602f-41c6-a5a4-0116879cd18c
newestrn = @reaction_network rn3 begin
            v, A + D --> 2D
           end v

# ╔═╡ 5c96f8dc-1932-4aea-a38b-9bd5b1a0e9d5
@named composed_network = compose(basern, [newrn, newestrn])

# ╔═╡ 900dd51e-4f5f-4f09-87c5-433433724fb0
md"Conceptually, we can think of `newrn` and `newestrn` as sub-systems in a tree, where `basern` is the root. We can visualize this like:"

# ╔═╡ 44eec599-dd11-4543-bf13-11d880edddf6
plot(TreePlot(composed_network), method=:tree, fontsize=12, nodeshape=:ellipse)

# ╔═╡ 570586d3-0192-4f2c-acbc-3cc55cb3c834
md"We can get the sub-systems of `composed_network` like"

# ╔═╡ e6e2572b-3730-43e5-bcff-1ed2ce6f9bc0
subsys = ModelingToolkit.get_systems(composed_network)

# ╔═╡ a09d9983-faca-47e0-b2e2-a4a29b435bcb
nameof(subsys[1]), nameof(subsys[2])

# ╔═╡ 79038b34-3fdb-4f56-bff6-cedf9518ef1f
md"
- When converting a `ReactionSystem` to ODEs, SDEs or jump processes, Catalyst first flattens these systems, merging all the systems together and renaming sub-system species and parameters to include their system name.
- This ensures they are treated as distinct species/parameters)
"

# ╔═╡ b6556575-1c48-42ab-b902-280c1239874d
flattened_rn = Catalyst.flatten(composed_network)

# ╔═╡ 4248e9d3-6415-4cf8-b20a-3e86e0b4e5f6
species(flattened_rn)

# ╔═╡ 6eeb56a3-da5c-4432-b939-a009138f2576
parameters(flattened_rn)

# ╔═╡ 8575c033-32fd-4c1b-88b8-6c372e449a87
reactions(flattened_rn)

# ╔═╡ ad102575-2957-4a20-9390-4650666fd335
md"The Catalyst [compositional modeling](https://catalyst.sciml.ai/dev/tutorials/compositional_modeling/) tutorial in the Catalyst docs has much more details on how to compose systems together and share variables between systems."

# ╔═╡ 137ef85f-2f1c-4890-9885-4a894e26dbfa
md"## _*The Repressilator via Compositional Modeling*_
The repressilator is a model of three genes where each gene represses the next in the sequence. 

Let's see how we can use compositional modeling to build and simulate the repressilator.
"

# ╔═╡ b30d2470-fd4a-4f6b-a7e6-9cdfb6ed3be1
function repressed_gene(; R, name)
	rn = @reaction_network $name begin
		hillr($R,α,K,n), ∅ --> m
        (δ,γ), m <--> ∅
        β, m --> m + P
        μ, P --> ∅
    end α K n δ γ β μ

	# here we set default values for the parameters in rn
	p = (:α => .5, :K => 40, :n => 2, :δ => log(2)/120,
         :γ => 5e-3, :β => log(2)/6, :μ => log(2)/60)

	# here we set default values for the initial conditions in rn
	setdefaults!(rn, p)
	u₀ = [:m => 0., :P => 0.0]
	setdefaults!(rn, u₀)

	rn
end

# ╔═╡ 4308a174-67c8-4356-85d3-c7c41e0df85b
md"
- We first declare one variable, correpsonding to the protein that will repress the first gene. 
- In the repressilator this would be the protein made by the third gene.
- If we assume we will call the third gene `G3` then the protein in the system made by `repressed_gene` will be named `G3₊P` in our final system.
  - When flattening subsystems modeling toolkit replaces `G3.P` for accessing the `P` state of the `G3` system with a state named `G3₊P` in the flattened model.
  "

# ╔═╡ ba5ddec7-96c5-4fe1-94cc-07de17f4ce6c
@variables G3₊P(t)

# ╔═╡ 32c21d52-6501-468f-9c31-2e2681497851
@named G1 = repressed_gene(; R=ParentScope(G3₊P))

# ╔═╡ 185fe648-c485-4ab9-bbfc-3e3d8496ad39
md"Here `ParentScope(G3₊P)` tells ModelingToolkit that the variable `G3₊P` is coming from one level higher in the system hierarchy."

# ╔═╡ bfe77918-a87c-48b7-8d09-9f4bf4182523
@named G2 = repressed_gene(; R=ParentScope(G1.P))

# ╔═╡ fbe9c85d-3e8c-4ab0-99a8-5271d00cb86e
@named G3 = repressed_gene(; R=ParentScope(G2.P))

# ╔═╡ 41d06379-00a2-468f-a951-6d4f432a841b
@named repressilator = ReactionSystem(t; systems=[G1,G2,G3])

# ╔═╡ ee7afe3c-4989-4dd5-8706-a9ee523a36b7
plot(TreePlot(repressilator), method=:tree, fontsize=12, nodeshape=:ellipse)

# ╔═╡ 095a5aff-4aa9-4138-bec5-0f8ee5d37399
md"Let's make the first gene have some protein initially, and then simulate the model:"

# ╔═╡ 9329d5b9-d9f3-4bbf-a6de-24cea133ccbb
let
	u0 = [G1.P => 20.0]
	tspan = (0.0, 10000.0)
	oprob = ODEProblem(repressilator, u0, tspan)
	sol = solve(oprob, Tsit5())
	plot(sol)
end

# ╔═╡ 88350aaf-530d-49ef-a2d2-c1ce2eeb226a
md"## _*Application: Hodgkin-Huxley Equation as a Constraint ODE*_
- In many applications one may want to have coupled ODEs or algebraic equations associated with a chemical reaction network. 
- These can be used in Catalyst via constraint systems.
- Constraint systems can be ModelingToolkit `ODESystem`s or `NonlinearSystem`s.

For example, let's build a simple Hodgkin-Huxley model for a single neuron, with the voltage, `V(t)`, included as a constraint `ODESystem`.

We first specify the transition rates for the three gating variables, ``m(t)``, ``n(t)``, and ``h(t)`` 
```math
s' \xleftrightarrow[\beta_s(V(t))]{\alpha_s(V(t))} s, \quad s \in \{m,n,h\}
```
where 
- ``m``, ``n``, and ``h``, are gating variables that determine the fraction of active (open) or inactive (``m' = 1 - m``, ``n' = 1 -n``, ``h' = 1 - h``) receptors.

The transition rate functions, which depend on the voltage, ``V(t)``, are then
"

# ╔═╡ 1b16a199-ffae-45e5-8528-21f05c4d3432
begin 
	function αₘ(V) 
		theta = (V + 45) / 10
		IfElse.ifelse(theta == 0.0, 1.0, theta/(1 - exp(-theta)))
	end
	βₘ(V) = 4*exp(-(V + 70)/18)
	
	αₕ(V) = .07 * exp(-(V + 70)/20)
	βₕ(V) = 1/(1 + exp(-(V + 40)/10))
	
	function αₙ(V)
		theta = (V + 60) / 10
		IfElse.ifelse(theta == 0.0, .1, .1*theta / (1 - exp(-theta)))
	end
	βₙ(V) = .125 * exp(-(V + 70)/80)
end

# ╔═╡ 10c8bca1-7f9f-457e-a19f-fc60d358317c
md"
- We now declare the symbolic variable, `V(t)`, that will represent voltage.
- We tell Catalyst not to generate an equation for it from the reactions we list, using the `isbcspecies` metadata. 
- This label tells Catalyst an ODE or nonlinear equation for `V(t)` will be provided in a constraint system.

Aside: `bcspecies` means a boundary condition species, a terminology from SBML.
"

# ╔═╡ 5390f023-b554-43be-b4f7-1738a3a2cab3
@variables V(t) [isbcspecies=true]

# ╔═╡ fea64be0-3d41-4714-8690-fbd5df42c246
hhrn = @reaction_network hhmodel begin
	(αₙ($V),βₙ($V)), n′ <--> n
	(αₘ($V),βₘ($V)), m′ <--> m
	(αₕ($V),βₕ($V)), h′ <--> h
end

# ╔═╡ b87af1c6-02b6-4981-9da4-c850c25d10eb
md"Next we create a `ModelingToolkit.ODESystem` to store the equation for `dV/dt`"

# ╔═╡ 0aa27dc1-8331-4d6a-8cf7-b15c91300458
voltageode = let
	@parameters C=1.0 ḡNa=120.0 ḡK=36.0 ḡL=.3 ENa=45.0 EK=-82.0 EL=-59.0 I₀=0.0
	@variables m(t) n(t) h(t)
	I = I₀* sin(2*pi*t/30)^2 

	Dₜ = Differential(t)
	eqs = [Dₜ(V) ~ -1/C * (ḡK*n^4*(V-EK) + ḡNa*m^3*h*(V-ENa) + ḡL*(V-EL)) + I/C]
	@named voltageode = ODESystem(eqs, t)
end

# ╔═╡ 07537b15-68e7-4642-8e8c-d5001b564101
md"Notice, we included an applied current, `I`, that we will use to perturb the system and create action potentials. For now we turn this off by setting its amplitude, `I₀`, to zero.


Finally, we couple this ODE into the reaction model as a constraint system:"

# ╔═╡ 2ff15451-de82-4547-abdd-c22fc498d352
@named hhmodel = extend(voltageode, hhrn)

# ╔═╡ f7cdaa82-9343-4b0e-a7b1-139109961951
md"
- `hhmodel` is now a `ReactionSystem` that is coupled to an internal constraint `ODESystem` that stores `dV/dt`. 
- Let's now solve to steady-state, as we can then use these resting values as an initial condition before applying a current to create an action potential."

# ╔═╡ 0c2af32c-70fc-4ed9-9673-601dc4191f24
hhsssol = let
	tspan = (0.0, 50.0)
	u₀ = [:V => -70, :m => 0.0, :h => 0.0, :n => 0.0, 
		  :m′ => 1.0, :n′ => 1.0, :h′ => 1.0]
	oprob = ODEProblem(hhmodel, u₀, tspan)
	sol = solve(oprob, Rosenbrock23())	
end

# ╔═╡ 7bcb71ea-5c83-4680-804e-79d2754d7b91
plot(hhsssol, vars=V)

# ╔═╡ b05cfc61-e815-4684-8014-bda19c9aef0e
u_ss = hhsssol.u[end]

# ╔═╡ e9232601-8b28-4e4c-abc0-dc5a8c46b7f4
md"Finally, starting from this resting state let's solve the system when the amplitude of the stimulus is non-zero and see if we get action potentials"

# ╔═╡ b16c5e12-32ac-4e40-8185-35efb1221b50
let
	tspan = (0.0, 50.0)
	@unpack I₀ = hhmodel
	oprob = ODEProblem(hhmodel, u_ss, tspan, [I₀ => 10.0])
	sol = solve(oprob)
	plot(sol, vars=V, legend=:outerright)
end

# ╔═╡ 7a065cca-ffac-47cb-bc29-c00237d6936a
md"## _*Appendix*_"

# ╔═╡ 34a4ecb6-c5ec-4b18-ad98-51769d5e8a0a
md"Let's set some default initial values for the voltage and gating variables."

# ╔═╡ 322c1baf-eaee-4342-a097-6206c95cb7c7
# let
# 	V₀ = -70
# 	setdefaults!(hhrn, [:noff => βₙ(V₀)/(αₙ(V₀)+βₙ(V₀)), :n => αₙ(V₀)/(αₙ(V₀)+βₙ(V₀)),
# 				    	:moff => βₘ(V₀)/(αₘ(V₀)+βₘ(V₀)), :m => αₘ(V₀)/(αₘ(V₀)+βₘ(V₀)),
# 						:hoff => βₕ(V₀)/(αₕ(V₀)+βₕ(V₀)), :h => αₕ(V₀)/(αₕ(V₀)+βₕ(V₀)),
# 						:V => V₀])
# end

# ╔═╡ a812e534-9df6-477d-9bf4-afba5aa4a45c


# ╔═╡ Cell order:
# ╟─192af62c-08fc-11ed-35a9-bb2a5d18c32f
# ╟─7cb53ff0-1751-4940-a4b3-bec68e81a005
# ╟─32897144-3033-45bb-a878-f72f6882df20
# ╟─146fee58-6588-4bfc-94d9-cfa7302b59a0
# ╠═3177b2e3-eefd-4faf-b292-2c7c6f0bbe22
# ╠═d35195ef-de7a-4f4c-8a38-ca9b4ec1b742
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
# ╠═bf439579-ab39-4e3e-9c1d-0d71ac7d8b9e
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
# ╟─b84c5bff-c4e3-41ed-8686-6253b4b87673
# ╠═bcdd2649-93c4-4d20-a4f6-e4437b8b842b
# ╟─c3be82d6-54f5-4876-8583-4f834417b0bf
# ╠═56c15ca4-1a26-4f7f-9d4a-dbe138614f03
# ╠═b4ad6737-55d8-417c-bc31-04898807089b
# ╟─7f6e982b-8710-4d1a-84ee-0ed77e568a8c
# ╠═306bcb0c-5134-4e50-a72a-5d559a1726d7
# ╠═a5277543-582d-4958-80f4-1f43857601a0
# ╠═63fd312f-50da-46e5-934a-6deacfbd474a
# ╟─cbb87922-da86-494e-8e46-7ecbad29f432
# ╠═ff6cd1cf-3753-4f93-94e1-dbbaccf97aa8
# ╟─1da19c48-9532-4e48-9bd5-ce968e3de493
# ╠═edd33cf5-b00e-4b13-bdfe-175cb01f6822
# ╟─92c7f19d-9299-488e-87e9-76c68fe26405
# ╠═e7a89226-4757-42c1-825a-b117ae439c1e
# ╟─9252661d-6fb8-47fb-805d-895f228df4e2
# ╠═03ce0b66-6860-4f8b-9f22-ec9ffaec60d1
# ╠═4430dc04-3a31-408c-85d9-1efe36c088d1
# ╠═0fe3d571-acf6-462c-a7f3-07e65ac0849f
# ╟─b34914dc-bdbb-4ba0-b992-b4295cf0c642
# ╠═60263535-9291-483f-b28e-fc8c26f4019f
# ╟─121b11c2-aadc-4e73-871a-d6c79c8fa3d9
# ╠═95e3bf53-d644-480a-9c75-602501657f7e
# ╟─15dede87-f769-4d16-ac76-95bf4d8777cd
# ╟─a4bd43d3-2b45-4a64-b17c-ca188d1d7d34
# ╠═a891e0ce-0cf5-4e5e-9ffb-b8cf257dd750
# ╟─82ac4773-4320-4683-85d1-8bdc0f5a77a3
# ╠═5535c9cc-bb87-4046-8621-93214829a012
# ╟─ed06ac37-98bb-4a28-81c7-29885405fa3b
# ╠═771e0080-a6d4-4de8-a7ec-3bb3949fc57f
# ╟─09d51be7-0972-4701-bfd7-471b207f518d
# ╟─b6dcccae-1aa3-4778-9b5b-bd0d8858df95
# ╠═fac7b265-4802-4fa2-87b8-a08a57e6ec6c
# ╠═becb9021-3e24-47a3-8893-3e5a55242958
# ╠═489b6939-5505-4bf7-8e5f-1a48ae4e6fe4
# ╟─09f90ebf-4fa8-4cca-9fb5-dcfc64084d69
# ╠═0d736c32-602f-41c6-a5a4-0116879cd18c
# ╠═5c96f8dc-1932-4aea-a38b-9bd5b1a0e9d5
# ╟─900dd51e-4f5f-4f09-87c5-433433724fb0
# ╠═44eec599-dd11-4543-bf13-11d880edddf6
# ╟─570586d3-0192-4f2c-acbc-3cc55cb3c834
# ╠═e6e2572b-3730-43e5-bcff-1ed2ce6f9bc0
# ╠═a09d9983-faca-47e0-b2e2-a4a29b435bcb
# ╟─79038b34-3fdb-4f56-bff6-cedf9518ef1f
# ╠═b6556575-1c48-42ab-b902-280c1239874d
# ╠═4248e9d3-6415-4cf8-b20a-3e86e0b4e5f6
# ╠═6eeb56a3-da5c-4432-b939-a009138f2576
# ╠═8575c033-32fd-4c1b-88b8-6c372e449a87
# ╟─ad102575-2957-4a20-9390-4650666fd335
# ╟─137ef85f-2f1c-4890-9885-4a894e26dbfa
# ╠═b30d2470-fd4a-4f6b-a7e6-9cdfb6ed3be1
# ╟─4308a174-67c8-4356-85d3-c7c41e0df85b
# ╠═ba5ddec7-96c5-4fe1-94cc-07de17f4ce6c
# ╠═32c21d52-6501-468f-9c31-2e2681497851
# ╟─185fe648-c485-4ab9-bbfc-3e3d8496ad39
# ╠═bfe77918-a87c-48b7-8d09-9f4bf4182523
# ╠═fbe9c85d-3e8c-4ab0-99a8-5271d00cb86e
# ╠═41d06379-00a2-468f-a951-6d4f432a841b
# ╠═ee7afe3c-4989-4dd5-8706-a9ee523a36b7
# ╟─095a5aff-4aa9-4138-bec5-0f8ee5d37399
# ╠═9329d5b9-d9f3-4bbf-a6de-24cea133ccbb
# ╟─88350aaf-530d-49ef-a2d2-c1ce2eeb226a
# ╠═1b16a199-ffae-45e5-8528-21f05c4d3432
# ╟─10c8bca1-7f9f-457e-a19f-fc60d358317c
# ╠═5390f023-b554-43be-b4f7-1738a3a2cab3
# ╠═fea64be0-3d41-4714-8690-fbd5df42c246
# ╟─b87af1c6-02b6-4981-9da4-c850c25d10eb
# ╠═0aa27dc1-8331-4d6a-8cf7-b15c91300458
# ╟─07537b15-68e7-4642-8e8c-d5001b564101
# ╠═2ff15451-de82-4547-abdd-c22fc498d352
# ╟─f7cdaa82-9343-4b0e-a7b1-139109961951
# ╠═0c2af32c-70fc-4ed9-9673-601dc4191f24
# ╠═7bcb71ea-5c83-4680-804e-79d2754d7b91
# ╠═b05cfc61-e815-4684-8014-bda19c9aef0e
# ╟─e9232601-8b28-4e4c-abc0-dc5a8c46b7f4
# ╠═b16c5e12-32ac-4e40-8185-35efb1221b50
# ╟─7a065cca-ffac-47cb-bc29-c00237d6936a
# ╟─34a4ecb6-c5ec-4b18-ad98-51769d5e8a0a
# ╠═322c1baf-eaee-4342-a097-6206c95cb7c7
# ╠═a812e534-9df6-477d-9bf4-afba5aa4a45c
