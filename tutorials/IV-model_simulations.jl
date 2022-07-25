### A Pluto.jl notebook ###
# v0.19.10

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
	sprob_ns = SDEProblem(pos_feedback_model_noise_scaling,u0,tspan,p_ns;noise_scaling=(@parameters η)[1]);
	ssol_ns = solve(sprob_ns,ImplicitEM());
	plot(ssol_ns)
end

# ╔═╡ 9dcfeb31-61ca-4a36-93bf-e20e9c852811
md"
## For high noise, the adaptive time stepper might fail.
Here, adaptive time-stepping can be turned off.
"

# ╔═╡ 84f026d3-927d-435d-a8fb-7984412b486a
let
	p_ns = [p; :η => 6.0]
	sprob_ns = SDEProblem(pos_feedback_model_noise_scaling,u0,tspan,p_ns;noise_scaling=(@parameters η)[1]);
	ssol_ns = solve(sprob_ns,ImplicitEM());
	plot(ssol_ns)
end

# ╔═╡ 8c2baa2d-67f7-49ec-9ce4-688b3bf88a08
md"
Here, \"adaptive=false\" turns of adaptiveness, and \"dt=0.001\" sets the time step size.
"

# ╔═╡ 9f536f9d-15b1-4499-be35-debf0048ed5c
let
	p_ns = [p; :η => 6.0]
	sprob_ns = SDEProblem(pos_feedback_model_noise_scaling,u0,tspan,p_ns;noise_scaling=(@parameters η)[1]);
	ssol_ns = solve(sprob_ns,ImplicitEM();adaptive=false,dt=0.001);
	plot(ssol_ns)
end

# ╔═╡ 2a5584f3-4a33-439f-a0de-6ddb767937cf
md"
For low dt, it might make sense to save the solution less frequently (by using the \"saveat=0.1\" option).
"

# ╔═╡ 1e37edbb-58e9-4f79-95ab-3dab20b9a586
let
	p_ns = [p; :η => 6.0]
	sprob_ns = SDEProblem(pos_feedback_model_noise_scaling,u0,tspan,p_ns;noise_scaling=(@parameters η)[1]);
	ssol_ns = solve(sprob_ns,ImplicitEM();adaptive=false,dt=0.001,saveat=0.1);
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

# ╔═╡ 2c51fc8d-083b-45a0-a663-26a269f3cc07
let
	p_activation = [:i => 0.0, :v0 => 1.0, :v => 10.0, :K => 100.0, :n => 2, :d => 0.1]
	sprob_activation = SDEProblem(pos_feedback_model,u0,tspan,deepcopy(p_activation))
	
	prob_func(prob,i,repeat) = prob
	eprob = EnsembleProblem(sprob_activation,prob_func=prob_func,safetycopy=false);
	
	esol = solve(eprob,ImplicitEM();trajectories=15,callback=activation_cb)
	plot(esol)
end

# ╔═╡ a2094021-cbfe-4d8f-8559-a37955c2eb34
md"
This looks good, but the callback changes the parameter value, which is then used in subsequent simulation (already before time 20). We need to fix this.
"

# ╔═╡ f6b373a5-888e-4fec-8872-b7dcdf5558b0
esol = let
	p_activation_en = [:i => 0.0, :v0 => 1.0, :v => 10.0, :K => 100.0, :n => 2, :d => 0.1]
	sprob_activation = SDEProblem(pos_feedback_model,u0,tspan,deepcopy(p_activation))
	
	prob_func(prob,i,repeat) = remake(prob,p=deepcopy(sprob_activation.p));
	eprob = EnsembleProblem(sprob_activation,prob_func=prob_func,safetycopy=false);
	solve(eprob,ImplicitEM();trajectories=15,callback=activation_cb)
end;

# ╔═╡ bad36227-f32d-4398-a5a1-62e92bf3105f
plot(esol)

# ╔═╡ 0275e54e-4e80-4015-8634-4ffb8f74785a
md"
The \"linealpha\" option can be useful when plotting a large number of solutions.
"

# ╔═╡ 4ffba508-491f-4065-a7b7-73dcd1f80ab5
plot(esol,linealpha=0.6)

# ╔═╡ 4c1e82c2-717c-45e1-8c0d-1ccaf0fa2f19
md"
A subset of the solutions can be plotted using the \"idxs\" argument:
"

# ╔═╡ f126a0a3-8c4d-4dee-9fec-b98f54aa860b
plot(esol,linealpha=0.6,idxs=1:5)

# ╔═╡ 0e10eabc-81dd-4187-9d6c-8491307a95e1
md"
## Further investigation of the solution object
"

# ╔═╡ a701bccd-f184-4281-b569-355f6ab9ee10
md"
We simulate a simple two-state model:
"

# ╔═╡ 573950d4-a321-420b-b425-d4a90c56185f
sol = let
	two_state_model = @reaction_network begin
		(k1,k2), X1 <--> X2
	end k1 k2
	sprob = SDEProblem(two_state_model,[:X1 => 5.0, :X2 => 5.0],(0.,10.0),[:k1 => 1.0, :k2 => 0.5])
	solve(sprob,ImplicitEM())
end;	

# ╔═╡ d0cd5b36-6849-45c9-8122-3f5a63ba22f7
plot(sol)

# ╔═╡ 5a418c56-6862-4dcd-87e3-c36a9477ab63
md"
One can select which variables to plot using the \"vars\" argument.
"

# ╔═╡ 28be06f5-b60b-41b8-a5ad-6deb69b9fd13
plot(sol,vars=[2])

# ╔═╡ d42dda26-f043-4ba0-9276-80e0b2bcf20c
md"
The solution can be accessed through the \"u\" field:
"

# ╔═╡ ccc64756-87a3-4911-ad99-2443f8828c81
sol.u

# ╔═╡ e9570862-4021-44f4-b8e4-a374489f1301
md"
Each element is a vector with the state at that time point.
"

# ╔═╡ 97ea0edd-ff92-4a83-bcbc-5ce5a3e90bb3
md"
And the time object through the \"t\" field.
"

# ╔═╡ f773f53e-6407-4775-b2aa-fca261d2bed6
sol.t

# ╔═╡ 6713871b-027c-40ad-876e-c54491c60fbf
md"
!!! note 
*Later we will show how one can easily use symbolic variables to access the solution and make plots, avoiding the need to know the integer index of a variable.*"

# ╔═╡ ece97e85-012d-4bc0-8933-3a27a05ede28
md"
## Jacobian for simulations
The jacobian can be automatically computed, given the correct argument.

For quick simulations, this is not a big concern, but can be for larger ones.
"

# ╔═╡ 9a531c2e-474c-425c-9757-ed3b77c2388a
begin
	brusselator_model = @reaction_network begin
	    A, ∅ → X
	    1, 2X + Y → 3X
	    B, X → Y
	    1, X → ∅
	end A B;
	u0_b = [:X => 5.0, :Y => 0.5]
	tspan_b = (0.,200000.0)
	p_b = [:A => 1.0, :B =>3.0]
end

# ╔═╡ 346a8fe9-75e0-448f-809d-ab54e63f681e
md"
Simulation not using the Jacobian:
"

# ╔═╡ ba19a7e7-8c99-477f-bdf0-f0c3cfcea605
let
	oprob = ODEProblem(brusselator_model,u0_b,tspan_b,p_b)
	@time solve(oprob)
end;

# ╔═╡ c880916d-f3e2-40c4-809a-cef081e93b88
md"
Simulation using the Jacobian:
"

# ╔═╡ 6ce9d36b-b934-43be-a978-3bf4e759a1bc
let
	oprob = ODEProblem(brusselator_model,u0_b,tspan_b,p_b;jac=true)
	@time solve(oprob)
end;

# ╔═╡ 29e081a9-0df3-4522-aff0-bf062486ba9a
md"
Simulation using the sparse Jacobian:
"

# ╔═╡ 48d084d5-09a6-4c67-b95f-3840889a3db6
let
	oprob = ODEProblem(brusselator_model,u0_b,tspan_b,p_b;jac=true,sparse=true)
	@time solve(oprob)
end;

# ╔═╡ 6ff919f2-39f7-4cbd-a6d2-6e57ca3fc3bc
md"
For a large model (BCR), using the Jacobian makes simulations twice as fast. Setting it to sparse provides a speed up of almost 2 orders of magnitude!
"

# ╔═╡ f146c1ba-7132-4c80-bc00-1b62f416ecad
md"
# Exercise: Simulate birth-death model.
Simulate the following birth-death model for the following parameter sets. What happens when fluctuations become large? What happens when concentrations are low (Hint: for low copy number simulations, a Gillespie approach is better than a CLE approach).

Next, for p1, p2, p3. Compare the effect of selecting these parameter sets to that of scaling the noise in the SDE problem.
"

# ╔═╡ 37ced29b-b58b-4e90-bf57-4ed9b39621a5
begin
	bd_model = @reaction_network begin
		(p,d), 0 <--> X
	end p d;
	p1 = [10.0, 1.0]
	p2 = 10*[10.0, 1.0]
	p3 = 100*[10.0, 1.0]
	p4 = [1.0, 1.0]
	p5 = 10*[1.0, 1.0]
end

# ╔═╡ Cell order:
# ╟─b67e5b7f-3664-4ae2-bf0c-aaadb0165dca
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
# ╟─9dcfeb31-61ca-4a36-93bf-e20e9c852811
# ╠═84f026d3-927d-435d-a8fb-7984412b486a
# ╟─8c2baa2d-67f7-49ec-9ce4-688b3bf88a08
# ╠═9f536f9d-15b1-4499-be35-debf0048ed5c
# ╟─2a5584f3-4a33-439f-a0de-6ddb767937cf
# ╠═1e37edbb-58e9-4f79-95ab-3dab20b9a586
# ╟─83ec221d-5ed0-4616-bf26-4e17af26cfe4
# ╠═c4757730-4cff-415c-991d-b6e917ac6288
# ╠═2235f660-7601-4e04-ab4d-77e33c09ff22
# ╟─70a154c3-a3d6-4082-9a53-e3869943af5e
# ╠═2c51fc8d-083b-45a0-a663-26a269f3cc07
# ╟─a2094021-cbfe-4d8f-8559-a37955c2eb34
# ╠═f6b373a5-888e-4fec-8872-b7dcdf5558b0
# ╠═bad36227-f32d-4398-a5a1-62e92bf3105f
# ╟─0275e54e-4e80-4015-8634-4ffb8f74785a
# ╠═4ffba508-491f-4065-a7b7-73dcd1f80ab5
# ╟─4c1e82c2-717c-45e1-8c0d-1ccaf0fa2f19
# ╠═f126a0a3-8c4d-4dee-9fec-b98f54aa860b
# ╟─0e10eabc-81dd-4187-9d6c-8491307a95e1
# ╟─a701bccd-f184-4281-b569-355f6ab9ee10
# ╠═573950d4-a321-420b-b425-d4a90c56185f
# ╠═d0cd5b36-6849-45c9-8122-3f5a63ba22f7
# ╟─5a418c56-6862-4dcd-87e3-c36a9477ab63
# ╠═28be06f5-b60b-41b8-a5ad-6deb69b9fd13
# ╟─d42dda26-f043-4ba0-9276-80e0b2bcf20c
# ╠═ccc64756-87a3-4911-ad99-2443f8828c81
# ╟─e9570862-4021-44f4-b8e4-a374489f1301
# ╟─97ea0edd-ff92-4a83-bcbc-5ce5a3e90bb3
# ╠═f773f53e-6407-4775-b2aa-fca261d2bed6
# ╟─6713871b-027c-40ad-876e-c54491c60fbf
# ╟─ece97e85-012d-4bc0-8933-3a27a05ede28
# ╠═9a531c2e-474c-425c-9757-ed3b77c2388a
# ╟─346a8fe9-75e0-448f-809d-ab54e63f681e
# ╠═ba19a7e7-8c99-477f-bdf0-f0c3cfcea605
# ╟─c880916d-f3e2-40c4-809a-cef081e93b88
# ╠═6ce9d36b-b934-43be-a978-3bf4e759a1bc
# ╟─29e081a9-0df3-4522-aff0-bf062486ba9a
# ╠═48d084d5-09a6-4c67-b95f-3840889a3db6
# ╟─6ff919f2-39f7-4cbd-a6d2-6e57ca3fc3bc
# ╟─f146c1ba-7132-4c80-bc00-1b62f416ecad
# ╠═37ced29b-b58b-4e90-bf57-4ed9b39621a5
