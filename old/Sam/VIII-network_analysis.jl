### A Pluto.jl notebook ###
# v0.19.10

using Markdown
using InteractiveUtils

# ╔═╡ ce709a57-d5d4-40de-916d-a5f032c54a2d
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, ModelingToolkit, DifferentialEquations, Plots, HomotopyContinuation
end

# ╔═╡ 402cfff0-0b8f-11ed-21dc-2dc56de9b6c9
html"<button onclick=present()>Present</button>"

# ╔═╡ 301e9dc2-273e-4d9d-a78f-48e4c1e458ff
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ 46ee8822-e319-404a-ba81-fd6bdd6ab373
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

# ╔═╡ 88cc932c-1f75-4a20-a2a4-aa7ae6a16562
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

# ╔═╡ 943caf53-af1b-4162-9bf1-99436a2214de
md"# Network Analysis
Catalyst offers a variety of functionality for analyzing and modifying chemical reaction networks, including:
- Detection and elimination of conservation laws.
- Calculation of stoichiometry matrices and rate law vectors.
- Calculation of reaction-complex network representations.
    - Calculation of linkage classes, network deficiency, reversibility, weak reversibility, etc...

Many of these are explored in the [Network Analysis in Catalyst](https://catalyst.sciml.ai/dev/tutorials/reaction_network_representation/) tutorial.

Let's see how we can calculate and use conservation laws via Catalyst.
"

# ╔═╡ f4bdb9f4-9a25-497e-bbc4-34833c0772f4
md"## _*Conservation Laws*_
Suppose we try to use HomotopyContinuation.jl to find the steady states of"

# ╔═╡ 22c586d0-e32c-4721-b8f8-5baadaf89c49
abtoc = @reaction_network begin
	(k₊,k₋), A + B <--> C
end k₊ k₋

# ╔═╡ e326ebdb-7641-4441-8cfa-f53df409b8d6
setdefaults!(abtoc, [:A => 100.0, :B => 80.0, :C => 20.0, :k₊ => 10.0, :k₋ => 1.0]);

# ╔═╡ 1647f9ac-fe88-454b-9a3f-97a3a00d3f52
ns = convert(NonlinearSystem, abtoc)

# ╔═╡ 2bda0427-1fe5-4736-ba09-eb348122f0b6
# for convenience
const MT = ModelingToolkit;

# ╔═╡ a8bba048-a1dd-44ab-a5a0-2f5f1cbc5b89
let	
	ps = MT.parameters(ns) .=> MT.varmap_to_vars([], MT.parameters(ns); defaults=MT.defaults(ns))
	subs = Dict(ps)
	new_eqs = map(eq -> substitute(eq.rhs, subs), equations(ns))
	sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))
end

# ╔═╡ 5aeee861-ef71-4e6d-99e0-83ee14736a7d
conservationlaws(abtoc)

# ╔═╡ 3f209388-8a58-4a12-a6ab-08cc17abb03f
string.(conservedequations(abtoc))

# ╔═╡ a41ae7a9-566b-4672-bb08-c9981e7809c4

string.(conservationlaw_constants(abtoc))

# ╔═╡ 63605f19-08f7-473a-b44f-3ff8a338f063
ns2 = convert(NonlinearSystem, abtoc; remove_conserved=true);

# ╔═╡ 1aedb584-d9af-4ef9-abf9-f3afe6eb597b
string(equations(ns2))

# ╔═╡ 742da7b2-1585-471b-bfd4-77df66776a71
ss = let	
	MT = ModelingToolkit

	# this gets the values for the reaction parameters and the conservation constants, evaluating the 
	# latter from the initial condition we specified
	ps = MT.parameters(ns2) .=> MT.varmap_to_vars([], MT.parameters(ns2); defaults=MT.defaults(ns2))
	
	subs = Dict(ps)
	new_eqs = map(eq -> substitute(eq.rhs, subs), equations(ns2))
	sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))
	vcat(filter(s -> s[1] >= 0, sols)...)
end

# ╔═╡ 668c8c04-02e0-4e9c-b676-7daa16e3c4c2
let
	oprob = ODEProblem(abtoc, [], (0.0, .03), []; remove_conserved=true)
	cons_sol = DifferentialEquations.solve(oprob, Tsit5())
	
	p = plot(cons_sol, legend=:outerright, size = (800,600))
	scatter!(p, cons_sol.t[end] * ones(length(ss)), ss, xlim=(0.0,1.1*cons_sol.t[end]), 
			 ylim=(0.0,100.0), label="steady-state")

	# note that we can still plot and manipulate the other variables too
	@unpack B,C = abtoc
	ps = plot(cons_sol, vars=[B,C], legend=:outerright)

	# let's calculate B+C which is conserved and plot it
	tv = range(0.0, .03, length=200)
	conserved_quant = [cons_sol(t, idxs=B) + cons_sol(t, idxs=C) for t in tv]
	pconslaw = plot(tv, conserved_quant, label="B(t) + C(t)", ylim=(90,110), legend=:outerright)
	
	plot(p, ps, pconslaw, layout=(3,1))
end

# ╔═╡ Cell order:
# ╟─402cfff0-0b8f-11ed-21dc-2dc56de9b6c9
# ╟─301e9dc2-273e-4d9d-a78f-48e4c1e458ff
# ╟─46ee8822-e319-404a-ba81-fd6bdd6ab373
# ╠═ce709a57-d5d4-40de-916d-a5f032c54a2d
# ╠═88cc932c-1f75-4a20-a2a4-aa7ae6a16562
# ╟─943caf53-af1b-4162-9bf1-99436a2214de
# ╟─f4bdb9f4-9a25-497e-bbc4-34833c0772f4
# ╠═22c586d0-e32c-4721-b8f8-5baadaf89c49
# ╠═e326ebdb-7641-4441-8cfa-f53df409b8d6
# ╠═1647f9ac-fe88-454b-9a3f-97a3a00d3f52
# ╠═2bda0427-1fe5-4736-ba09-eb348122f0b6
# ╠═a8bba048-a1dd-44ab-a5a0-2f5f1cbc5b89
# ╠═5aeee861-ef71-4e6d-99e0-83ee14736a7d
# ╠═3f209388-8a58-4a12-a6ab-08cc17abb03f
# ╠═a41ae7a9-566b-4672-bb08-c9981e7809c4
# ╠═63605f19-08f7-473a-b44f-3ff8a338f063
# ╠═1aedb584-d9af-4ef9-abf9-f3afe6eb597b
# ╠═742da7b2-1585-471b-bfd4-77df66776a71
# ╠═668c8c04-02e0-4e9c-b676-7daa16e3c4c2
