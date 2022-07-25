### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 5bac6466-4673-4d24-94dc-e43232942149
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, DifferentialEquations, Plots
end

# ╔═╡ 0794bd0a-80cd-476d-89dc-5cf9a8b3c289
html"<button onclick=present()>Present</button>"

# ╔═╡ 1641fd62-09c8-11ed-1891-7b343f6e889a
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

# ╔═╡ e557e806-3cfc-41fd-9bc8-c59851880a9b
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ f0389c62-1830-4875-920b-0d7567432750
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

# ╔═╡ c49d8402-0e67-41b6-b93d-a0e24a3902d1
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

# ╔═╡ 04b67692-4431-4508-bb79-31aa82065c8c
md"
# Introduction to basic Catalyst usage
Here we briefly introduce how to create and simulate a Catalyst model, before more thoroughly introducing the package.

A simple birth-death model (a single species being both produced and degraded) will be used for this example).
"

# ╔═╡ 454b973d-448f-41ee-a2ad-07acc5bdd669
md"
#### Creating the model
"

# ╔═╡ 8b5d701d-40ed-4699-b4fc-43521123d4ea
bd_model = @reaction_network begin
	b, 0 --> X
	d, X --> 0
end b d

# ╔═╡ fcdc330d-a027-40cc-9349-e8062891754c
md"
## Making an ODE simulation
"

# ╔═╡ 242235c2-e77a-45a1-ad53-95d948cd30d8
begin
	u0 = [:X => 1.0];
	tspan = (0.0,50.0)
	p = [:b => 1.0, :d => 0.1]
	oprob = ODEProblem(bd_model,u0,tspan,p)
end;

# ╔═╡ 790ed9f4-84ff-4284-a4a0-21c5b76cb4ea
begin
	osol = solve(oprob)
	plot(osol)
end

# ╔═╡ f692fc4a-34c5-449b-818b-95f8a5ff9d07
md"
## Making an SDE simulation
"

# ╔═╡ 74f2b081-bded-4632-bb4f-265c78f8ff85
begin
	sprob = SDEProblem(bd_model,u0,tspan,p);
	ssol = solve(sprob)
	plot(ssol)
end

# ╔═╡ ea3a847b-32ab-46fd-b02c-8dff06e3eb66
md"
## Making a Gillespie simulation
"

# ╔═╡ 6d0e3d7f-a028-4f21-acbc-d22ccdcc29af
begin
	u0_gill = [:X => 1];
	dprob = DiscreteProblem(bd_model,u0_gill,tspan,p)
	jprob = JumpProblem(bd_model,dprob,Direct());
end;

# ╔═╡ 5c66f199-1c9e-41c8-93cb-f49c04e827c5
begin
	gsol = solve(jprob,SSAStepper())
	plot(gsol)
end

# ╔═╡ Cell order:
# ╟─0794bd0a-80cd-476d-89dc-5cf9a8b3c289
# ╟─1641fd62-09c8-11ed-1891-7b343f6e889a
# ╟─e557e806-3cfc-41fd-9bc8-c59851880a9b
# ╟─f0389c62-1830-4875-920b-0d7567432750
# ╟─5bac6466-4673-4d24-94dc-e43232942149
# ╟─c49d8402-0e67-41b6-b93d-a0e24a3902d1
# ╟─04b67692-4431-4508-bb79-31aa82065c8c
# ╟─454b973d-448f-41ee-a2ad-07acc5bdd669
# ╠═8b5d701d-40ed-4699-b4fc-43521123d4ea
# ╟─fcdc330d-a027-40cc-9349-e8062891754c
# ╠═242235c2-e77a-45a1-ad53-95d948cd30d8
# ╠═790ed9f4-84ff-4284-a4a0-21c5b76cb4ea
# ╟─f692fc4a-34c5-449b-818b-95f8a5ff9d07
# ╠═74f2b081-bded-4632-bb4f-265c78f8ff85
# ╟─ea3a847b-32ab-46fd-b02c-8dff06e3eb66
# ╠═6d0e3d7f-a028-4f21-acbc-d22ccdcc29af
# ╠═5c66f199-1c9e-41c8-93cb-f49c04e827c5
