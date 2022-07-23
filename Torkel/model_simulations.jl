### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 956875d4-4b94-47cc-91a2-1830c08fa09f
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, DifferentialEquations, Plots
end

# ╔═╡ b67e5b7f-3664-4ae2-bf0c-aaadb0165dca
html"<button onclick=present()>Present</button>"

# ╔═╡ 1ad036f8-0a7e-11ed-1f0c-a9dfa382bf47
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

# ╔═╡ 0de20f4f-5123-4924-b037-196ebe474f8a
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ fe4ad7fa-b363-424f-b074-49c8a834280a
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

# ╔═╡ 1ba42536-51bc-46b8-ad2d-b2be28f6966d
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

# ╔═╡ 4d5a20d1-aeea-42da-b924-01fcb06fafb0
md"
# Simulation of Catalyst models
"

# ╔═╡ ab14af3d-35cb-4e9d-9fa0-6bdcd557fab0
md"
#### Creating the model
A positive feedback loop where the single species (*X*) activates its own production.
"

# ╔═╡ 5a6bfba4-d8c8-4d1d-9ea4-db1971361f31
pos_feedback_model = @reaction_network begin
    v0 + hill(i*X,v,K,n), ∅ → X
    d, X → ∅
end i v0 v K n d

# ╔═╡ 9427c782-63a0-41ec-be07-e26cd5e65228
md"
## Making an ODE simulation
"

# ╔═╡ 9a10ab17-cea5-4d68-88da-2d5d0568c5c6
begin
	u0 = [:X => 5.0]
	tspan = (0.0,100.0)
	p = [:i => 1.0, :v0 => 1.0, :v => 10.0, :K => 100.0, :n => 2, :d => 0.1]
	oprob = ODEProblem(pos_feedback_model,u0,tspan,p)
end;

# ╔═╡ 5dbd36aa-5659-4e2c-b642-c8782fbff67d
begin
	osol = solve(oprob,Rosenbrock23())
	plot(osol)
end

# ╔═╡ 7ed8c738-2441-4198-a4ce-abe96bdf8b2a
plot(osol,xguide="Time",yguide="Concentration",framestyle=:box,label="",lw=4)

# ╔═╡ 93abfa3f-b178-433c-9a9f-3febe20b7457
md"
## Making an SDE simulation
"

# ╔═╡ 370b4310-618e-45ac-9c03-75e142ddcccd
begin
	sprob = SDEProblem(pos_feedback_model,u0,tspan,p);
	ssol = solve(sprob,ImplicitEM())
	plot(ssol)
end

# ╔═╡ d866c9ed-c8e8-4ec3-b009-87e556907358
md"
## Noise scaling in the CLE interpretation
Under certain circumstances, it might be desirable to tune the amount of noise used in the CLE SDE interpretation. We can add a parameter that linearly tunes the noise amplitude:
"

# ╔═╡ 44ff366b-adb3-427d-a71f-22e88b5ef02d
pos_feedback_model_noise_scaling = @reaction_network begin
    v0 + hill(i*X,v,K,n), ∅ → X
    d, X → ∅
end i v0 v K n d η;

# ╔═╡ 78135a37-8ddd-4281-a1d7-2fea75286383
begin
	p_ns = [p; :η => 0.1]
	oprob_ns = SDEProblem(pos_feedback_model_noise_scaling,u0,tspan,p_ns;noise_scaling=(@parameters η)[1]);
	ssol_ns = solve(oprob_ns,ImplicitEM());
	plot(ssol_ns)
end

# ╔═╡ 83ec221d-5ed0-4616-bf26-4e17af26cfe4
md"
## Callbacks permit mid-simulation events
"

# ╔═╡ c4757730-4cff-415c-991d-b6e917ac6288
begin
	affect!(integrator) = integrator.p[1]=3.0
	activation_cb = PresetTimeCallback([20.0],affect!)
end;

# ╔═╡ 2235f660-7601-4e04-ab4d-77e33c09ff22
begin
	p_activation = [:i => 0.0, :v0 => 1.0, :v => 10.0, :K => 100.0, :n => 2, :d => 0.1]
	sprob_activation = SDEProblem(pos_feedback_model,u0,tspan,deepcopy(p_activation))
	ssol_activation = solve(sprob_activation,ImplicitEM(),callback=activation_cb,adaptive=false,dt=0.0001)
	plot(ssol_activation)
end

# ╔═╡ 70a154c3-a3d6-4082-9a53-e3869943af5e
md"
## Ensemble Problems permit Monte-Carlo simulations
"

# ╔═╡ f6b373a5-888e-4fec-8872-b7dcdf5558b0
begin
	sprob_activation_ens = SDEProblem(pos_feedback_model,u0,tspan,deepcopy(p_activation))
	prob_func(prob,i,repeat) = remake(prob,p=deepcopy(sprob_activation_ens.p));
	eprob = EnsembleProblem(sprob_activation_ens,prob_func=prob_func,safetycopy=false);
	esol = solve(eprob,ImplicitEM();trajectories=15,callback=activation_cb)
	plot(esol)
end

# ╔═╡ 4ffba508-491f-4065-a7b7-73dcd1f80ab5
plot(esol,linealpha=0.6)

# ╔═╡ ac0598c9-43bb-4499-8151-ec76263970cd
solve(eprob,ImplicitEM(),EnsembleThreads();trajectories=5,callback=activation_cb);

# ╔═╡ Cell order:
# ╟─b67e5b7f-3664-4ae2-bf0c-aaadb0165dca
# ╟─1ad036f8-0a7e-11ed-1f0c-a9dfa382bf47
# ╟─0de20f4f-5123-4924-b037-196ebe474f8a
# ╟─fe4ad7fa-b363-424f-b074-49c8a834280a
# ╟─956875d4-4b94-47cc-91a2-1830c08fa09f
# ╟─1ba42536-51bc-46b8-ad2d-b2be28f6966d
# ╟─4d5a20d1-aeea-42da-b924-01fcb06fafb0
# ╟─ab14af3d-35cb-4e9d-9fa0-6bdcd557fab0
# ╠═5a6bfba4-d8c8-4d1d-9ea4-db1971361f31
# ╟─9427c782-63a0-41ec-be07-e26cd5e65228
# ╠═9a10ab17-cea5-4d68-88da-2d5d0568c5c6
# ╠═5dbd36aa-5659-4e2c-b642-c8782fbff67d
# ╠═7ed8c738-2441-4198-a4ce-abe96bdf8b2a
# ╟─93abfa3f-b178-433c-9a9f-3febe20b7457
# ╠═370b4310-618e-45ac-9c03-75e142ddcccd
# ╟─d866c9ed-c8e8-4ec3-b009-87e556907358
# ╠═44ff366b-adb3-427d-a71f-22e88b5ef02d
# ╠═78135a37-8ddd-4281-a1d7-2fea75286383
# ╟─83ec221d-5ed0-4616-bf26-4e17af26cfe4
# ╠═c4757730-4cff-415c-991d-b6e917ac6288
# ╠═2235f660-7601-4e04-ab4d-77e33c09ff22
# ╟─70a154c3-a3d6-4082-9a53-e3869943af5e
# ╠═f6b373a5-888e-4fec-8872-b7dcdf5558b0
# ╠═4ffba508-491f-4065-a7b7-73dcd1f80ab5
# ╠═ac0598c9-43bb-4499-8151-ec76263970cd
