### A Pluto.jl notebook ###
# v0.19.10

using Markdown
using InteractiveUtils

# ╔═╡ 74931442-a55b-4331-b570-9c3a352e5cf4
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, ModelingToolkit, DifferentialEquations, Plots, StochasticDiffEq, JumpProcesses, OrdinaryDiffEq
end

# ╔═╡ 47323c44-09fe-11ed-0e8c-1b3df7e3b610
html"<button onclick=present()>Present</button>"

# ╔═╡ 5eb4e150-9594-4de4-9016-259b4d84e694
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ 5b4f1bcf-cf28-4531-a273-5869fe7bc82b
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

# ╔═╡ f5188c47-4c58-46a6-8e93-2b9b7b5d575a
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

# ╔═╡ 1e81937b-b71d-48a7-bdce-1a3f612c35be
md"# What's going on under the hood?
Let's look at what happens when we create various problem types, like `ODEProblem`, from a Catalyst reaction model. We'll focus on
- Default reaction rate laws and `ModelingtToolkit.ODESystem`
- Alternative rate laws via the `combinatoric_ratelaws=false` option.
- Other ModelingToolkit system types.
"

# ╔═╡ 25f947a8-0352-4588-b7e6-ff715c367d1e
md"## _*What happens when we create an `ODEProblem` from a reaction network?*_
Suppose we have a reaction system with
- ``L`` species, ``(S_1,\dots,S_L)``.
Consider the reaction
```math
m_1 S_1 + m_2 S_2 + \dots m_{L} S_{L} \overset{k}{\to} n_1 S_1 + n_2 S_2 + \dots n_L S_L
```
We call
- ``(m_1,\dots,m_L)`` the *substrate* stoichiometric coefficients.
- ``(n_1,\dots,n_L)`` the *product* stoichiometric coefficients.
- ``(\nu_1,\dots,\nu_L) = (n_1-m_1,\dots,n_L-m_L)`` the *net* stoichiometric coefficients.
- ``\nu_{\ell}`` gives the change in ``S_{\ell}`` when the reaction occurs.

An `ODEProblem` corresponds to a Catalyst-generated reaction rate equation (RRE) ODE model for a given reaction system. For the reaction above we would have 
```math
\begin{align*}
\frac{dS_{\ell}}{dt} &= \nu_{\ell} R(S_1,\dots,S_L) \\
&= \nu_{\ell} k \prod_{i = 1}^L \frac{(S_{i})^{m_i}}{m_i!}
\end{align*}
```
Here ``R(S_1,\dots,S_L)`` is called the *rate law* of the reaction.

In Catalyst we could create a reaction network with this reaction like
```julia
rn = @reaction_network begin
	k, m_1*S_1 + ... + m_L*S_L --> n_1 S_1 + ... + n_L S_L
end k
```
for `m_1`, ..., `m_L` and `n_1`, ..., `n_L` user-specified positive integers. 
"

# ╔═╡ 489608bc-e9da-487d-a65b-3ff4ef24b2f4
md"## _*Can we get at the ODEs that Catalyst generates?*_
Yes! Just as Catalyst gives a symbolic representation of a reaction network, in creating an `ODEProblem` we actually (internally) first create a `ModelingToolkit.ODESystem`, which is a symbolic model for a system of ODEs.

- ModelingToolkit has many system types, corresponding to symbolic mathematical models. 
- Each is an instance of a `ModelingToolkit.AbstractSystem`, with a unified set of core accessor functions (see ModelingToolkit.jl).
- This include `ODESystem`, `NonlinearSystem`, `SteadyStateSystem`, `SDESystem`, and `JumpSystem`.


An example of an ODESystem:
"

# ╔═╡ 0d85b799-20b3-47d3-a760-6f44e9aca2e8
rn = @reaction_network begin
	k, 2*X + 3*Y --> 4Z
end k

# ╔═╡ 0f975934-cfd2-4272-9f52-5f57815540f5
rn isa ReactionSystem

# ╔═╡ dedcca44-ebbf-4b99-ba66-dc2fe37dfca5
ReactionSystem <: ModelingToolkit.AbstractSystem

# ╔═╡ ed82e9e0-ac63-4771-9e99-52c687d94de4
md"So we see that all the `@reaction_network` macro is doing is creating a `ReactionSystem`, which is a type of `ModelingToolkit.AbstractSystem`."

# ╔═╡ 78b92ef8-d622-486a-b9f6-e9ff57c2ce0e
md"Previously we created an `ODEProblem` like"

# ╔═╡ 2abf6f02-8e81-42d8-8439-dbbb5c733d9a
u₀ = [:X => 1.0, :Y => 1.0, :Z => 1.0];

# ╔═╡ 76a43a31-5a02-47c7-af5b-53e46b972df7
p = [:k => 1.0];

# ╔═╡ f33f84c9-5830-4afc-b620-f2fef6ee7ccb
tspan = (0.0, 10.0);

# ╔═╡ 818314bd-9d8d-4f08-b4d0-9521cede2674
oprob = ODEProblem(rn, u₀, tspan, p);

# ╔═╡ 5d96198d-cd71-4cff-9db2-cf536448ac2c
md"Under the hood, Catalyst calling `ODEProblem(rn, u₀, tspan, p)` first creates an `ODESystem` with the symbolic ODEs, and then calls ModelingToolkit's `ODEProblem` on this system, i.e. what Catalyst is doing internally is"

# ╔═╡ 21c1982e-0799-4b34-ab0b-940ed8952586
osys = convert(ODESystem, rn)

# ╔═╡ 447663db-a95f-47fd-a6f9-e5c4e7318325
md"The symbolic ODEs one would actually type in creating an `ODESystem` directly are:"

# ╔═╡ 9606facf-3e43-4c58-bdb7-ca283308cd2b
string.(equations(osys))

# ╔═╡ 09e3e79d-cd0e-429e-bff1-bc43696877eb
md"We'll see more about creating ODEs when we look at coupling them to `ReactionSystem`s later in the workshop."

# ╔═╡ 74dd38e4-9eb8-4c46-83bc-75538c2ab66d
md"Here we see the combinatorial factors that arise from the factorials within each ODE!

ModelingToolkit doesn't support the `Symbol`-based method by which we specified `u₀` and `p`, so that we can't use them in `ODEProblem(osys,...)`. However, we can convert them to the form ModelingToolkit uses and then pass them in as follows"

# ╔═╡ a64ab3f1-8e14-44b1-b126-e68dbd9e5a5c
u₀MT = symmap_to_varmap(osys, u₀);

# ╔═╡ f9c99fba-0cd8-4047-ba7a-27e4c08a0b0d
pMT = symmap_to_varmap(osys, p);

# ╔═╡ 8b737401-7d9a-4d83-a627-8b43ff95dac6
oprob2 = ODEProblem(osys, u₀MT, tspan, pMT);

# ╔═╡ cc4e327c-ecc3-4639-b9e1-cca2a8efde31
md"`oprob` and `oprob2` ultimately encode the same set of ODEs. 

Let's check this by solving and plotting them:"

# ╔═╡ abd7d4db-0be2-436f-9532-50891c14c0ce
let
	sol = solve(oprob, Tsit5())
	sol2 = solve(oprob2, Tsit5())
	plot(plot(sol, title="oprob", legend=:right), plot(sol2, title="oprob2", legend=:right))
end

# ╔═╡ 96047157-4d4c-429d-82a1-679e0d9872e6
md"## _*Aside: what if one does not want the factorial scaling?*_
In many cases the rate constant, ``k``, already contains the factorial terms that appear in the denominator of the rate law, i.e. one wants the ODEs

```math
\begin{align*}
\frac{dS_{\ell}}{dt} = \nu_{\ell} k \prod_{i = 1}^L (S_{i})^{m_i}.
\end{align*}
```
This can be controlled via the `combinatoric_ratelaws` keyword argument
"

# ╔═╡ 98cd78c0-126c-42d1-80ef-a54e4fc151e6
oprob_nocr = ODEProblem(rn, u₀, tspan, p; combinatoric_ratelaws = false);

# ╔═╡ c0bc648e-0499-4225-a397-cf2513bd742a
md"or"

# ╔═╡ c0f88f5d-16e6-47c4-bfe3-2d5c7fdf6ce2
osys_nocr = convert(ODESystem, rn; combinatoric_ratelaws = false)

# ╔═╡ 4ec5b213-ecd9-4866-af68-51c3fe7ab23e
md"Compare to the previous case:"

# ╔═╡ 1787ba78-fd36-4dd0-85ca-b29c64adf858
osys

# ╔═╡ 2f898587-ea3f-40c2-bb8c-e03fc4f8f83e
md"Let's create the associated `ODEProblem` and show `oprob_nocr` and it give the same answer."

# ╔═╡ 1f2332be-7bec-4ae9-a99f-d0e6b5206bc1
oprob2_nocr = ODEProblem(osys_nocr, u₀MT, tspan, pMT);

# ╔═╡ 0506a612-0cff-429e-9a80-5ff7b770637f
let
	sol = solve(oprob_nocr, Tsit5())
	sol2 = solve(oprob2_nocr, Tsit5())
	plot(plot(sol, title="oprob_nocr", legend=:right), plot(sol2, title="oprob2_nocr", legend=:right))
end

# ╔═╡ 57588fcd-a0c5-4be9-9687-56bdd15a6d35
md"So we see we either 
- Set `combinatoric_ratelaws=false` when creating an `ODEProblem` directly from a `ReactionSystem`.
- Set `combinatoric_ratelaws=false` when creating an `ODESystem` from a `ReactionSystem`.
to drop the factorial scaling terms.
"

# ╔═╡ 3596e588-9421-46f0-91ea-84e74b624172
md"## _*ODEs vs. SDEs vs. Jump Processes via Catalyst*_ 
For the abstract ``m_1 S_1 + \dots m_L S_L \overset{k}{\to} n_1 S_1 \dots n_L S_L`` reaction, Catalyst
generates the following models (with `combinatoric_ratelaws=true`, the default).

Reaction rate equation (RRE) ODEs, `ModelingToolkit.ODESystem`:
```math
\begin{align*}
\frac{dS_{\ell}}{dt} &= \nu_{\ell} R(S_1,\dots,S_L) = \nu_{\ell} k \prod_{i = 1}^L \frac{(S_{i})^{m_i}}{m_i!}
\end{align*}
```
Chemical Langevin Equation (CLE) SDEs, `ModelingToolkit.SDESystem`:
```math
\begin{align*}
dS_{\ell}(t) &= \nu_{\ell} R(S_1,\dots,S_L) dt + \nu_{\ell} \sqrt{R(S_1,\dots,S_L)} dW(t),
\end{align*}
```
where ``dW(t)`` is a standard Wiener Process (Brownian Motion) with mean zero and variance of one.

Stochastic Chemical Kinetics jump process, `ModelingToolkit.JumpSystem`:
```math
S_{\ell}(t) = S_{\ell}(0) + \nu_{\ell} Y\paren{\int_0^t a(S_1(s^-),\dots,S_L(s^-)) \, ds},
```
where ``Y(t)`` denotes a unit Poission counting process, and the jump process transition rate (i.e. propensity or rate law) is
```math
a(S_1(t),\dots,S_L(t)) = k \prod_{i=1}^L \frac{S_i (S_i-1) \dots (S_i - m_i + 1)}{m_i!}
```

So in general
"

# ╔═╡ c8ce323a-86b4-43ec-ba7f-3433c4f6ffc2
combinatoric_ratelaws = true;  # or false!

# ╔═╡ dbfe7922-f4a9-4030-bf44-01d23a7991e6
sprob = SDEProblem(rn, u₀, tspan, p; combinatoric_ratelaws);

# ╔═╡ 6e3d1343-c37c-4f2b-9faf-9f4f988c52e9
md"is equivalent to"

# ╔═╡ 0c7b2e3a-98b7-49e2-8c76-7825c2faaedc
sdesys = convert(SDESystem, rn; combinatoric_ratelaws);

# ╔═╡ 9373bc36-f8b1-4cf2-99e7-0a51988db24d
sprob2 = SDEProblem(sdesys, u₀MT, tspan, pMT);

# ╔═╡ 91acc7bb-6606-4f16-8511-a70e754e8b47
md"and"

# ╔═╡ 65d8c886-2cc6-481c-b033-9a9a890164ff
dprob = DiscreteProblem(rn, u₀, tspan, p);

# ╔═╡ 19ccc574-0f36-496f-abc4-c35821c4e662
jprob = JumpProblem(rn, dprob, Direct(); combinatoric_ratelaws);

# ╔═╡ 4fceb9bf-51fa-423a-992e-1966abbbad8e
md"is equivalent to "

# ╔═╡ 82bb6f0c-2a9a-4942-9ebc-d971eea13fb4
jsys = convert(JumpSystem, rn; combinatoric_ratelaws);

# ╔═╡ ef362bd6-18e8-42e1-a7a1-00e069f8c9c8
dprob2 = DiscreteProblem(jsys, u₀MT, tspan, pMT);

# ╔═╡ 5285aead-f558-4eee-a1f2-6a521ec4716d
jprob2 = JumpProblem(jsys, dprob2, Direct());

# ╔═╡ 48678636-1c58-46b3-a91c-8ce0aa1443db
md"#### See the [first Catalyst tutorial](https://catalyst.sciml.ai/dev/tutorials/using_catalyst/#Reaction-rate-laws-used-in-simulations) for more information!"

# ╔═╡ Cell order:
# ╟─47323c44-09fe-11ed-0e8c-1b3df7e3b610
# ╟─5eb4e150-9594-4de4-9016-259b4d84e694
# ╟─5b4f1bcf-cf28-4531-a273-5869fe7bc82b
# ╠═74931442-a55b-4331-b570-9c3a352e5cf4
# ╠═f5188c47-4c58-46a6-8e93-2b9b7b5d575a
# ╟─1e81937b-b71d-48a7-bdce-1a3f612c35be
# ╟─25f947a8-0352-4588-b7e6-ff715c367d1e
# ╟─489608bc-e9da-487d-a65b-3ff4ef24b2f4
# ╠═0d85b799-20b3-47d3-a760-6f44e9aca2e8
# ╠═0f975934-cfd2-4272-9f52-5f57815540f5
# ╠═dedcca44-ebbf-4b99-ba66-dc2fe37dfca5
# ╟─ed82e9e0-ac63-4771-9e99-52c687d94de4
# ╟─78b92ef8-d622-486a-b9f6-e9ff57c2ce0e
# ╠═2abf6f02-8e81-42d8-8439-dbbb5c733d9a
# ╠═76a43a31-5a02-47c7-af5b-53e46b972df7
# ╠═f33f84c9-5830-4afc-b620-f2fef6ee7ccb
# ╠═818314bd-9d8d-4f08-b4d0-9521cede2674
# ╟─5d96198d-cd71-4cff-9db2-cf536448ac2c
# ╠═21c1982e-0799-4b34-ab0b-940ed8952586
# ╟─447663db-a95f-47fd-a6f9-e5c4e7318325
# ╟─9606facf-3e43-4c58-bdb7-ca283308cd2b
# ╟─09e3e79d-cd0e-429e-bff1-bc43696877eb
# ╟─74dd38e4-9eb8-4c46-83bc-75538c2ab66d
# ╠═a64ab3f1-8e14-44b1-b126-e68dbd9e5a5c
# ╠═f9c99fba-0cd8-4047-ba7a-27e4c08a0b0d
# ╠═8b737401-7d9a-4d83-a627-8b43ff95dac6
# ╟─cc4e327c-ecc3-4639-b9e1-cca2a8efde31
# ╠═abd7d4db-0be2-436f-9532-50891c14c0ce
# ╟─96047157-4d4c-429d-82a1-679e0d9872e6
# ╠═98cd78c0-126c-42d1-80ef-a54e4fc151e6
# ╟─c0bc648e-0499-4225-a397-cf2513bd742a
# ╠═c0f88f5d-16e6-47c4-bfe3-2d5c7fdf6ce2
# ╟─4ec5b213-ecd9-4866-af68-51c3fe7ab23e
# ╠═1787ba78-fd36-4dd0-85ca-b29c64adf858
# ╟─2f898587-ea3f-40c2-bb8c-e03fc4f8f83e
# ╠═1f2332be-7bec-4ae9-a99f-d0e6b5206bc1
# ╠═0506a612-0cff-429e-9a80-5ff7b770637f
# ╟─57588fcd-a0c5-4be9-9687-56bdd15a6d35
# ╟─3596e588-9421-46f0-91ea-84e74b624172
# ╠═c8ce323a-86b4-43ec-ba7f-3433c4f6ffc2
# ╠═dbfe7922-f4a9-4030-bf44-01d23a7991e6
# ╟─6e3d1343-c37c-4f2b-9faf-9f4f988c52e9
# ╠═0c7b2e3a-98b7-49e2-8c76-7825c2faaedc
# ╠═9373bc36-f8b1-4cf2-99e7-0a51988db24d
# ╟─91acc7bb-6606-4f16-8511-a70e754e8b47
# ╠═65d8c886-2cc6-481c-b033-9a9a890164ff
# ╠═19ccc574-0f36-496f-abc4-c35821c4e662
# ╟─4fceb9bf-51fa-423a-992e-1966abbbad8e
# ╠═82bb6f0c-2a9a-4942-9ebc-d971eea13fb4
# ╠═ef362bd6-18e8-42e1-a7a1-00e069f8c9c8
# ╠═5285aead-f558-4eee-a1f2-6a521ec4716d
# ╟─48678636-1c58-46b3-a91c-8ce0aa1443db
