### A Pluto.jl notebook ###
# v0.19.10

using Markdown
using InteractiveUtils

# ╔═╡ 89855cb1-eaac-4d35-b36a-9d317a4e6e1c
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, ModelingToolkit, DifferentialEquations, Plots, PlutoUI, DiffEqCallbacks
end

# ╔═╡ dada26bd-9564-4957-ba26-0c351ad59ecc
begin 
	using DiffEqProblemLibrary.JumpProblemLibrary
	JumpProblemLibrary.importjumpproblems()
end;

# ╔═╡ ceeace90-0b74-11ed-028d-316c8b46be70
html"<button onclick=present()>Present</button>"

# ╔═╡ 2522dd53-a78e-4fdf-a65e-8667c03d0054
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ 6dff89b6-2da5-44d8-91d5-48eaf7c4538b
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

# ╔═╡ 54a6ebb2-dcd3-4b2a-a376-e2eaa3dd7775
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

# ╔═╡ 698fcefd-bd81-4132-9707-aef8088e552f
md"# Jump Process Simulations
Stochastic simulation algorithm (SSA) simulations of stochastic chemical kinetics jump process models can be carried out by JumpProcesses.jl, a component in DifferentialEquations.jl. 

We'll now see some considerations when trying to run such simulations including:
1. Performance tuning:    
    - reducing memory use
    - different solvers
2. Using callbacks.

As an illustrating example, let's work with a toggle switch model for two genes that repress each other
"

# ╔═╡ fd3fb9d2-8960-48c9-be43-4a05d7fdaf5f
toggleswitch = @reaction_network toggleswitch begin
    hillr(P₂,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    μ, P₁ --> ∅
    μ, P₂ --> ∅
end α K n δ γ β μ

# ╔═╡ 2498ecb8-e8a0-46cd-bdea-86e60745ac10
begin 
	p = [:α => .15, :K => 40, :n => 2, :δ => log(2)/120, :γ => 5e-3, :β => 20*log(2)/120, :μ => log(2)/60]	
	u₀ = [:m₁ => 0, :m₂ => 0, :P₁ => 0, :P₂ => 0]
	tspan = (0., 100.)

	dprob = DiscreteProblem(toggleswitch, u₀, tspan, p)
	jprob = JumpProblem(toggleswitch, dprob, Direct())   
	sol = solve(jprob, SSAStepper()) 
	plot(sol, legend=:outerright)
end

# ╔═╡ d81f817c-ac77-4a3b-84fd-71f0eed41b3f
md"Note that here
1. `Direct()` indicates we want to use Gillespie's Direct Method as our SSA.
    - One can choose a number of alternative *aggregators* (in the JumpProcesses.jl language) as the SSA for selecting the next reaction time and type.
2. As we have a problem with no explicit time-dependencies, we can use the `SSAStepper()` time-stepper, which is optimized for such problems.
"

# ╔═╡ b22171d7-2173-4acf-a204-c2388e6ea09f
md"## _*How can we reduce memory use?*_
Notice, we aren't really seeing the state switching of the toggle switch above -- we need to run our simulation for a longer time. 
- Naively doing this will lead to a very large amount of data getting saved as the default is to save the state before and after *every* jump!
    - One benefit to this is that indexing the solution object gives the *exact* solution at a given time.
- We can modulate this by using the `save_positions` kwarg to `JumpProblem` with the `saveat` kwarg to solve:
"

# ╔═╡ 75d98830-cf97-495a-a364-cc4c275e4d80
let
	
	tspan = (0.0, 100000.0)
	dprob = DiscreteProblem(toggleswitch, u₀, tspan, p)

	# turn off saving the state before and after every jump	
	jprob = JumpProblem(toggleswitch, dprob, Direct(); save_positions=(false,false))

	# save the solution every 10.0 seconds only
	sol = solve(jprob, SSAStepper(), saveat=10.0)
	
	plot(sol, legend=:outerright)
end

# ╔═╡ 67f08ff1-ec51-4df3-b351-c02fbc4e91e8
md"## _*What are `ModelingToolkit.JumpSystem`s?*_
Let's take a look at how jumps are represented in a `ModelingToolkit.JumpSystem`:
"

# ╔═╡ d8811849-ebd3-4c3e-a1ed-e9cf49d63a76
jsys = convert(JumpSystem, toggleswitch);

# ╔═╡ 596c8f91-3327-4f77-990e-397015b5dca9
typeof.(equations(jsys))

# ╔═╡ 01d4181f-0869-4eed-8a05-a9464aa9c3d9
md"So we see the jumps that make up `jsys` are a mix of `MassActionJump`s and `ConstantRateJump`s. What are these? 

Broadly speaking (see [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) for more details), from least to most performant we have
- `VariableRateJump`s, the most general type of jump.
- `ConstantRateJump`s, jumps where the rate function (i.e. propensity, intensity or transition rate function), is **constant** between jumps of the system.
    - They can not have an explicit time-dependency in the rate function, or involve state variables that change from continuous dynamics (including states that are changed via a `VariableRateJump`).
- `MassActionJump`s are a special, performance-optimized, type of `ConstantRateJump` that efficiently represents mass action reactions via the rate laws used in Catalyst.


Today we will focus only on models for which the first two apply.
- In this case the `JumpProcesses.SSAStepper()` time-stepper is an efficient method for advancing the system from jump to jump.

Luckily Catalyst handles these classifications for us!
"

# ╔═╡ 7820e2d0-71a4-4d3f-8262-a1a328375323
md"## _*What types of SSAs are available to use?*_
A number of recent, high-performance SSAs for simulating stochastic chemical kinetics models are available to use with Catalyst models.
- Catlayst automatically handles generating all needed input for these methods for you! 
- This includes creating reaction-reaction, reaction-species, and species-reactions dependency graphs.
- Recall, SSAs are called aggregators in JumpProcesses.jl.


The current SSAs available for mixed `ConstantRateJump` and `MassActionJump` problems are
- *Gillespie's direct method (`Direct`)*
- Gillespie's first reaction method (`FRM`)
- Gibson and Bruck's next reaction method (`NRM`)
- *Sorting direct method (`SortingDirect`)*
- Reject direct method (`RDirect`)
- *Composition-rejection direct method (`DirectCR`)*
- *Rejection-SSA (RSSA)*
- *Rejection-SSA with composition-rejection (RSSACR)*
"

# ╔═╡ 559bf785-2b63-44f0-b0e3-29a9eeb983a0
md"## _*What are some benchmarks for a variety of system types?*_
Let's start with a small autoregulatory model of a gene with a protein product that can feedback and turn off the gene. 

We load our model from the `DiffEqProblemLibrary:`
"

# ╔═╡ d29e235b-41f9-4817-be2e-d9815664af03
negexpress = prob_jump_dnarepressor.network

# ╔═╡ ca752ac8-e0d1-475a-a0bf-a2586d07d805
md"We previously found the median timing relative to Direct is:
"

# ╔═╡ 405eaa72-d197-4bdc-9b2f-a8a97f9a6ada
LocalResource("negexpress.svg")

# ╔═╡ 06580418-a38b-49bd-a718-dcb6a44368cf
md"Next we consider a more complicated multistate model from Gupta and Mendes, *An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical Systems*, Computation, 6 (9), 2018.

This model is representative of a protein that can undergo reactions across a variety of phosphorylation states:"

# ╔═╡ 89776745-2b5b-417d-9d6a-3a50934cba8b
multistate = prob_jump_multistate.network

# ╔═╡ ebf890db-1d54-46f6-9338-7d5a3aebf2d8
md"We find"

# ╔═╡ b16d7769-c3ec-4af2-8e1f-e4849f028197
LocalResource("multistate.svg")

# ╔═╡ df28af85-954c-48e6-b0f4-331a68f52412
md"Finally, we consider a much larger model of a B cell receptor (BCR) signaling pathway.
- Over 1000 species and 24,000 reactions.
- We also include the BioNetGen simulator in our results, which uses a C\+\+ based version of the `SortingDirect` method.

This model can be loaded into Catalyst via [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl) package, see the examples folder there.

We find
"

# ╔═╡ f34a0386-9824-41ce-9a86-87f2620df72d
LocalResource("bcr.png")

# ╔═╡ 6b8866f3-d7ec-4d1b-899f-61c4996d52e1
md"*More benchmarking results will be in our forthcoming Catalyst.jl paper!*"

# ╔═╡ b0af5cd8-ab73-409d-9086-b19c8a0194c7
md"## _*How can we create events in SSA simulations?*_
`SSAStepper` supports discrete events in the DifferentialEquations.jl notation.
- These are events that occur based on a condition being true, which is checked between each jump of the system.
- The can be used to change the state, change parameters, or to terminate a simulation.

Let's see a few example of using them!

Let's first look at the preceding gene network, and calculate the first time that the gene becomes inactivated:
"


# ╔═╡ e9f08ed5-2807-493a-b659-ea0a038dc6c4
negexpress

# ╔═╡ 39943d98-4711-4a37-aa7b-2823dbcc5595
firsttimeprob, firsttimecb = let
	# set some initial conditions for our system
	setdefaults!(negexpress, [:DNA => 1, :DNAR => 0, :P => 0, :mRNA => 0])

	# convert to a JumpSystem
	jsys = convert(JumpSystem, negexpress)
	
	# symbolic variables we'll use
	@unpack DNA, DNAR = jsys

	# get its integer index within the solution vector, u
	DNARidx = findfirst(isequal(DNAR), species(jsys))

	# the condition function to stop a simulation
	function first_DNAR_cond(u, t, integrator; DNARidx=DNARidx)
		u[DNARidx] == 1
	end

	# an affect! to stop the simulation when the condition is true
	function stop_sim(integrator)
		savevalues!(integrator, true)
		terminate!(integrator)
		nothing
	end

	# create the callback
	cb = DiscreteCallback(first_DNAR_cond, stop_sim)
	
	dprob = DiscreteProblem(jsys, [], (0.0, 1000.0), prob_jump_dnarepressor.rates)
	jprob = JumpProblem(jsys, dprob, Direct())
	jprob,cb
end

# ╔═╡ fb48eabe-ea6f-4746-8ec5-ab55111a8ead
first_DNAR_sol = solve(firsttimeprob, SSAStepper(); save_end = false, callback = firsttimecb);

# ╔═╡ f23a97db-4000-4e5f-b2d4-3182c399cb2d
plot(first_DNAR_sol)

# ╔═╡ bdca5464-9bac-475f-8b9e-8e9c17e6d5bf
md"Let's look at the distribution of such exit times:"

# ╔═╡ 44db3690-450f-4fee-a99f-54ad8a352a1f
function getexittimes(jprob, cb, N)
	times = zeros(N)
	for i in eachindex(times)
		sol = solve(jprob, SSAStepper(); save_end = false, callback = cb)
		times[i] = sol.t[end]
	end
	times
end

# ╔═╡ a6cd74e7-f1a6-4f8e-9514-bdeb7e0fc7b0
times = getexittimes(firsttimeprob, firsttimecb, 100000)

# ╔═╡ 6f4eb0ad-7253-4ed4-8781-004f77c31ed4
histogram(times; normalize=true, bins=100, legend=nothing, xlabel="exit time", ylabel="PDF(exit time)")

# ╔═╡ cd1cbbd8-3ab5-4827-8ea8-4b2112a962c6
md"## _*What about if we want to change a state or parameters?*_
This is easy to do too, however, we need to make sure to tell the solver that its internal state may need to be reset due to our event changing a species or state!

We can do this using the `reset_aggregated_jumps!` command within our `affect!` function.

As a simple example, let's suppose we have cancerous cells that are produced at a constant rate, can die, and we apply a treatment that reduces their production rate by an additional 50% each time it is applied.
"

# ╔═╡ 692d54d0-fba1-41c0-8a6c-87098f354a83
let 
	birthdeath = @reaction_network begin
		b, ∅ --> N
		d, N --> ∅
	end b d

	jsys = convert(JumpSystem, birthdeath)
	
	# get the parameter's integer index in jsys
	@unpack b = jsys
	bidx = findfirst(isequal(b), parameters(jsys))

	# set up initial condition and parameters for jsys
	p = symmap_to_varmap(jsys, [:b => 2000.0, :d => 5.0])
	u₀ = symmap_to_varmap(jsys, [:N => 20])

	# we'll use a PresetTimeCallback to apply five treatments
	treatment_times = collect(range(4.0, 12.0, step=2.0))
	function treatment_affect!(integrator)
		integrator.p[bidx] *= .5
		reset_aggregated_jumps!(integrator)
		nothing
	end
	cb = PresetTimeCallback(treatment_times, treatment_affect!)

	dprob = DiscreteProblem(jsys, u₀, (0.0, 20.0), p)
	jprob = JumpProblem(jsys, dprob, Direct())
	sol = solve(jprob, SSAStepper())
	p1 = plot(sol, title="no treatment")
	
	sol = solve(jprob, SSAStepper(); callback = cb)
	p2 = plot(sol, title="with treatment")
	scatter!(p2, treatment_times, sol(treatment_times; idxs=jsys.N), label="treatment times")
	plot(p1,p2, layout=(2,1))
end

# ╔═╡ Cell order:
# ╟─ceeace90-0b74-11ed-028d-316c8b46be70
# ╟─2522dd53-a78e-4fdf-a65e-8667c03d0054
# ╟─6dff89b6-2da5-44d8-91d5-48eaf7c4538b
# ╠═89855cb1-eaac-4d35-b36a-9d317a4e6e1c
# ╠═54a6ebb2-dcd3-4b2a-a376-e2eaa3dd7775
# ╟─698fcefd-bd81-4132-9707-aef8088e552f
# ╠═fd3fb9d2-8960-48c9-be43-4a05d7fdaf5f
# ╠═2498ecb8-e8a0-46cd-bdea-86e60745ac10
# ╟─d81f817c-ac77-4a3b-84fd-71f0eed41b3f
# ╟─b22171d7-2173-4acf-a204-c2388e6ea09f
# ╠═75d98830-cf97-495a-a364-cc4c275e4d80
# ╟─67f08ff1-ec51-4df3-b351-c02fbc4e91e8
# ╠═d8811849-ebd3-4c3e-a1ed-e9cf49d63a76
# ╠═596c8f91-3327-4f77-990e-397015b5dca9
# ╟─01d4181f-0869-4eed-8a05-a9464aa9c3d9
# ╟─7820e2d0-71a4-4d3f-8262-a1a328375323
# ╟─559bf785-2b63-44f0-b0e3-29a9eeb983a0
# ╠═dada26bd-9564-4957-ba26-0c351ad59ecc
# ╠═d29e235b-41f9-4817-be2e-d9815664af03
# ╟─ca752ac8-e0d1-475a-a0bf-a2586d07d805
# ╟─405eaa72-d197-4bdc-9b2f-a8a97f9a6ada
# ╟─06580418-a38b-49bd-a718-dcb6a44368cf
# ╠═89776745-2b5b-417d-9d6a-3a50934cba8b
# ╟─ebf890db-1d54-46f6-9338-7d5a3aebf2d8
# ╟─b16d7769-c3ec-4af2-8e1f-e4849f028197
# ╟─df28af85-954c-48e6-b0f4-331a68f52412
# ╟─f34a0386-9824-41ce-9a86-87f2620df72d
# ╟─6b8866f3-d7ec-4d1b-899f-61c4996d52e1
# ╟─b0af5cd8-ab73-409d-9086-b19c8a0194c7
# ╠═e9f08ed5-2807-493a-b659-ea0a038dc6c4
# ╠═39943d98-4711-4a37-aa7b-2823dbcc5595
# ╠═fb48eabe-ea6f-4746-8ec5-ab55111a8ead
# ╠═f23a97db-4000-4e5f-b2d4-3182c399cb2d
# ╠═bdca5464-9bac-475f-8b9e-8e9c17e6d5bf
# ╠═44db3690-450f-4fee-a99f-54ad8a352a1f
# ╠═a6cd74e7-f1a6-4f8e-9514-bdeb7e0fc7b0
# ╠═6f4eb0ad-7253-4ed4-8781-004f77c31ed4
# ╟─cd1cbbd8-3ab5-4827-8ea8-4b2112a962c6
# ╠═692d54d0-fba1-41c0-8a6c-87098f354a83
