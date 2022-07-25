### A Pluto.jl notebook ###
# v0.19.10

using Markdown
using InteractiveUtils

# ╔═╡ e3b90e09-1613-4e48-8860-b50ac1c47fca
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, BifurcationKit, Setfield, LinearAlgebra, HomotopyContinuation, Plots
end

# ╔═╡ 0dd82ade-f03f-4fa6-be3a-93be2aa84fef
html"<button onclick=present()>Present</button>"

# ╔═╡ 4eb379f0-1492-4c5f-8e9d-2704eedf9602
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ 0649a89d-77f1-42f9-96e7-7f61cca954ea
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

# ╔═╡ d214fb01-ce3a-44d5-a036-1d88c13a4e58
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

# ╔═╡ ab7d1177-f46d-4cc1-be60-9fa78631189d
md"
# Plotting bifurcation diagrams using BifurcationKit.jl
BifurcationKit allows us to compute (and plot) *bifurcation diagrams*. These show how the system's steady states are affected by the value of some parameter.
"

# ╔═╡ 8777ec41-71c3-43e9-993f-4e4069890dcf
md"
We will use the positive feedback model, which we already suspect might be bistable.
"

# ╔═╡ b064289b-14ae-4441-b871-369e1e41df8a
pos_feedback_model = @reaction_network begin
    v0 + hill(i*X,v,K,n), ∅ → X
    d, X → ∅
end i v0 v K n d;

# ╔═╡ b5483297-2580-43c5-b3bc-6e001ad938b0
md"
Next, we extract the ODE function in a form that BifurcationKit can appreciate.
"

# ╔═╡ 0fbfcbd4-fdb0-4b4a-bf5f-9254aa08638a
begin
	odefun = ODEFunction(convert(ODESystem,pos_feedback_model),jac=true)
	F = (u,p) -> odefun(u,p,0)      
	J = (u,p) -> odefun.jac(u,p,0)
	jet = BifurcationKit.getJet(F, J; matrixfree=false)
end;

# ╔═╡ f02d05b5-ec3f-4087-b4ff-fb9dc3c4d830
md"
We have to designate the parameter set for which we want to compute our bifurcation diagram. In addition, we need to select our bifurcation parameter, and over which range we wish to vary it.
"

# ╔═╡ 6df3870e-38fc-40a6-967a-02f0d37f443f
begin
	params = [0.5,0.1,10.0,0.5,2,1.0]
	bif_par_idx = 1
	p_span = (0.0,1.0)
end;

# ╔═╡ 1ed9f8a1-4661-4ed1-96b2-8d50ce68f06f
md"
We select the variable which we wish to plot in our bifurcation diagram (in this case there is only one alternative).
"

# ╔═╡ 27669090-bf88-4e1f-90ff-f3ea17277b46
plot_var_idx = 1;

# ╔═╡ 6064a994-a536-443f-b989-e9b9bace42cd
md"
To ensure the diagram is plotted across the whole parameter range, we need to modify our parameter set so that the target parameter is at the beginning of the interval.
"

# ╔═╡ 6a43f65e-5817-436a-8381-1dde21e606ae
begin\
	bif_params = copy(params)
	bif_params[bif_par_idx] = p_span[1]
end;

# ╔═╡ 03c706e8-6c92-43f1-b516-6672ae9ec041
md"
Finally, we need to provide a guess for the steady state(s) at the initial point of the bifurcation diagram. This could be achieved by running an ODE simulation. However, as exact values are not needed, we will simply make a rough guess.
"

# ╔═╡ f1244183-4f89-4281-90b7-3c9cdb6added
u0 = [0.1];

# ╔═╡ 0e8e295c-8800-49dc-8360-a5f57756d2fb
md"
Next, we need to set our continuation parameters. These are used when the bifurcation algorithm is computed.
"

# ╔═╡ 2282c4f1-50a1-402b-920d-723fdeffee18
opts_br = ContinuationPar(
		pMin=p_span[1], pMax=p_span[2],
		dsmax = 1e-2, dsmin = 1e-5, ds=1e-3, maxSteps=1000,
		detectBifurcation=3);

# ╔═╡ 9de0ecc6-9bc6-40b3-9978-48c8143d3a3a
md"
We can now compute the bifurcation diagram:
"

# ╔═╡ 822ab1af-b7f9-4204-a2aa-ab5da9e47e8e
bif_dia = bifurcationdiagram(jet..., u0, bif_params, (@lens _[bif_par_idx]), 2,
                (x,p,level)->setproperties(opts_br);
                tangentAlgo = BorderedPred(),
                recordFromSolution=(x, p) -> x[plot_var_idx], verbosity = 0, plot=false);

# ╔═╡ b394fa10-9cd5-4bcf-9791-0ee431fc70ec
md"
The bifurcation diagram shows how the concentration of X depend on the value of the parameter i. The thick lines denote stable steady states, and the thin one unstable ones. Bifurcation points (points where the quantity or quality of the steady states changes) are also displayed.
"

# ╔═╡ b532a718-df9c-4935-b23c-a21b7e9dfead
plot(bif_dia)

# ╔═╡ 105e49ae-dd43-4431-8006-226b9a1ab12a
md"
BifurcationKit is a powerful package, that is capable of a lot more. In addition, for some systems, additional considerations are required. Please consult https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/ for more details
"

# ╔═╡ bfaa59d1-95db-427f-aa2c-7fb0af0a0449
md"
!!! note 
This approach only works if the bifurcation diagram can be tracked as one continious line. If we reduce the parameter range (so that the first bifurcation point is not reached), the bistable region is not properly found. There are alternative methods for handling these cases.
"

# ╔═╡ 586f9f1e-1e53-4091-839f-b301e0683437
let
	p_span = (0.0,0.2)
	opts_br = ContinuationPar(
			pMin=p_span[1], pMax=p_span[2],
			dsmax = 1e-2, dsmin = 1e-5, ds=1e-3, maxSteps=1000,
			detectBifurcation=3);
	bif_dia = bifurcationdiagram(jet..., u0, bif_params, (@lens _[bif_par_idx]), 2,
	                (x,p,level)->setproperties(opts_br);
	                tangentAlgo = BorderedPred(),
	                recordFromSolution=(x, p) -> x[plot_var_idx], verbosity = 0, plot=false);
	plot(bif_dia)
end

# ╔═╡ 8b1b0583-85ea-404f-8703-4c0cc1e77662
md"
# Finding steady states using HomotopyContinuation.jl
This is a powerful package that can find all solutions to n multivariate polynomials of n variables. This is a large class of systems that includes most CRNs. 
"

# ╔═╡ 968f458f-5022-47cc-a667-038cd80bcdc8
md"
We use a model from the Wilhem (2009) paper (which demonstrates bistability in a small CRN). We declare the model and the parameter set for which we want to find the steady states.
"

# ╔═╡ ff2eaf0a-074a-4401-8a55-9d1756c29474
begin
	wilhelm_2009_model = @reaction_network begin
	    k1, Y --> 2X
	    k2, 2X --> X + Y
	    k3, X + Y --> Y
	    k4, X --> 0
	end k1 k2 k3 k4
	p = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]
end;

# ╔═╡ d3a3784c-76e9-4752-87a3-f02409ae339e
md"
Next, we need to extract the system which we need to solve to find the steady states. We also need to insert the parameter values. In the future, a helper function for this will likely be created for this and the next step.
"

# ╔═╡ c8b94894-0ca7-4997-866b-ea092ae29105
begin
	ns = convert(NonlinearSystem,wilhelm_2009_model)
	subs = Dict(Pair.(wilhelm_2009_model.ps.value,last.(p)))
	new_eqs = map(eq -> ModelingToolkit.unwrap(substitute(eq.rhs,subs)), ns.eqs.value)
end

# ╔═╡ 40058a11-195e-4443-a946-6be9dbabd7e6
md"
Finally, we use the \"as_polynomial\" function to read the equations as a polynomial, and solve them using homotopy continuation. In addition, we filter away solutions where any value is imaginary.
"

# ╔═╡ 939c4ca8-db48-4a08-bfc7-0632a0315a58
sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))

# ╔═╡ 11ed4d07-0104-4994-9fea-14eef6393aac
md"
While it is not the case for this CRN, we note that some solutions with negative species concentrations may still appear, in which case these needs also be filtered out.
"

# ╔═╡ c195dc55-4db9-44ed-ac8d-1fef9810811c
md"
# Excersise: Bifurcation Diagram for Wilhelm 2009 Model
Consider the following, slightly modified model from the previous section. Try computing the bifurcation diagram for the parameter k2, over the interval (1.0,15.0). Attempt to use homotopy Continuation to confirm it finds the same fixed points for various selections of k2. Hint: You might need to modify the \"maxSteps\" argument.
"

# ╔═╡ 580f6825-2f27-4187-a5c3-087479941bc4
begin
	wilhelm_2009_moddified_model = @reaction_network begin
		0.1, 0 --> X
	    k1, Y --> 2X
	    k2, 2X --> X + Y
	    k3, X + Y --> Y
	    k4, X --> 0
	end k1 k2 k3 k4
end;

# ╔═╡ 92789d15-12a1-4045-8ca0-c77d1e7f457b


# ╔═╡ Cell order:
# ╟─0dd82ade-f03f-4fa6-be3a-93be2aa84fef
# ╟─4eb379f0-1492-4c5f-8e9d-2704eedf9602
# ╟─0649a89d-77f1-42f9-96e7-7f61cca954ea
# ╟─e3b90e09-1613-4e48-8860-b50ac1c47fca
# ╟─d214fb01-ce3a-44d5-a036-1d88c13a4e58
# ╟─ab7d1177-f46d-4cc1-be60-9fa78631189d
# ╟─8777ec41-71c3-43e9-993f-4e4069890dcf
# ╠═b064289b-14ae-4441-b871-369e1e41df8a
# ╟─b5483297-2580-43c5-b3bc-6e001ad938b0
# ╠═0fbfcbd4-fdb0-4b4a-bf5f-9254aa08638a
# ╟─f02d05b5-ec3f-4087-b4ff-fb9dc3c4d830
# ╠═6df3870e-38fc-40a6-967a-02f0d37f443f
# ╟─1ed9f8a1-4661-4ed1-96b2-8d50ce68f06f
# ╠═27669090-bf88-4e1f-90ff-f3ea17277b46
# ╟─6064a994-a536-443f-b989-e9b9bace42cd
# ╠═6a43f65e-5817-436a-8381-1dde21e606ae
# ╟─03c706e8-6c92-43f1-b516-6672ae9ec041
# ╠═f1244183-4f89-4281-90b7-3c9cdb6added
# ╟─0e8e295c-8800-49dc-8360-a5f57756d2fb
# ╠═2282c4f1-50a1-402b-920d-723fdeffee18
# ╟─9de0ecc6-9bc6-40b3-9978-48c8143d3a3a
# ╠═822ab1af-b7f9-4204-a2aa-ab5da9e47e8e
# ╟─b394fa10-9cd5-4bcf-9791-0ee431fc70ec
# ╠═b532a718-df9c-4935-b23c-a21b7e9dfead
# ╟─105e49ae-dd43-4431-8006-226b9a1ab12a
# ╟─bfaa59d1-95db-427f-aa2c-7fb0af0a0449
# ╠═586f9f1e-1e53-4091-839f-b301e0683437
# ╟─8b1b0583-85ea-404f-8703-4c0cc1e77662
# ╟─968f458f-5022-47cc-a667-038cd80bcdc8
# ╠═ff2eaf0a-074a-4401-8a55-9d1756c29474
# ╟─d3a3784c-76e9-4752-87a3-f02409ae339e
# ╠═c8b94894-0ca7-4997-866b-ea092ae29105
# ╟─40058a11-195e-4443-a946-6be9dbabd7e6
# ╠═939c4ca8-db48-4a08-bfc7-0632a0315a58
# ╟─11ed4d07-0104-4994-9fea-14eef6393aac
# ╟─c195dc55-4db9-44ed-ac8d-1fef9810811c
# ╠═580f6825-2f27-4187-a5c3-087479941bc4
# ╟─92789d15-12a1-4045-8ca0-c77d1e7f457b
