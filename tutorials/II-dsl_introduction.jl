### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 073c9ca2-7504-42ec-ad9f-e9e3014303f3
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Catalyst, DifferentialEquations, Latexify, Plots
end

# ╔═╡ 37a68fe9-328d-4e10-84da-113dae3b1a9e
html"<button onclick=present()>Present</button>"

# ╔═╡ 3283a80a-09cc-11ed-2959-6525dd49b3fd
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

# ╔═╡ 6ec927cc-a67b-44ba-b750-6947a1518dbc
html"""<style>
main {
    max-width: 900px;
}
"""

# ╔═╡ 76af30dd-f995-4a48-b922-76c2e84b1693
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

# ╔═╡ 35be8093-beb7-45f9-bfa1-51e91e27292a
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

# ╔═╡ 0879683c-b797-4ba1-8e86-8a16006bbdc5
show_odes(rn) = latexify(convert(ODESystem,rn));

# ╔═╡ c78dbeaf-bda0-4be7-816b-553bf7300829
md"
# Creating models using the catalyst DSL
Catalyst provides a DSL (domain-specific language) to enable simple declaration of CRN models.
"

# ╔═╡ 344dbafb-0d22-40f4-a250-e16ac366a41a
md"
## Declaring a birth-death model
"

# ╔═╡ 66e4e1c9-3678-4235-8bc6-71c24b1cd1d5
md"
Here, each reaction event of the model is listed. First, the rate is declared, followed by the reaction (substrates and products separated by an arrow). All parameters are listed at the end.
"

# ╔═╡ 32ec973b-4e14-461a-bf30-7262265c0fe0
bd_model = @reaction_network begin
	b, 0 --> X
	d, X --> 0
end b d

# ╔═╡ 06d766b4-261b-46ac-a948-ea51ffe42572
md"
Latexify can be used to investigate the corresponding ODE model:
"

# ╔═╡ a03aec3b-82c6-4dc6-b046-20f7f7b04c5f
show_odes(bd_model)     # (Special function used for this presentation, not exported)

# ╔═╡ fb29bcf8-e30c-4ff5-bed5-82314d371204
md"
## Species stoichiometries
"

# ╔═╡ 620f072c-1d4b-4ae3-895e-c792471502fe
md"
A number before a species denotes its stoichiometry in a reaction.
"

# ╔═╡ c0ab7ba6-0437-44cc-85d9-af6d2753339a
dimerisation_model = @reaction_network begin
	kB, 2X --> X2
	kD, X2 --> 2X
end kB kD

# ╔═╡ 83b8f3fd-bff1-4588-b951-87863ff885b5
md"
Higher order reactions (stoichiometries > 1) are accounted for in the ODE.
"

# ╔═╡ 024163ae-0ad4-41fd-a48b-0ed85876d960
show_odes(dimerisation_model)   

# ╔═╡ 2feb2fe7-2c07-4d7e-9b19-e8e7d68f1e24
md"
Another (artificial) example:
"

# ╔═╡ b838469c-7658-41a5-af1c-272fe1ceff68
stoich_example_model = @reaction_network begin
	r, 2X + 3Y --> Z
end r;

# ╔═╡ fbdcc69a-6849-48c4-8e8c-de3c31bb77c2
show_odes(stoich_example_model)

# ╔═╡ 3f1c0b3c-f988-4180-b637-ba6e11f8ada8
md"
# Bundling of similar reactions
Many CRNs contain nearly identical reactions, these can often be bundles to provide more concise notation.
"

# ╔═╡ 7cfadb90-4283-4471-b707-373b01aa4fd8
md"
## Bi-directional reactions
Reversible reactions (typically with different rates in each direction) can be declared using a single line.
"

# ╔═╡ ad17c467-0852-4414-a2c4-df8b36329846
md"
The birth-death model essentially consists of a single, reversible, reaction.
"

# ╔═╡ a4b6ad01-2713-4927-89c2-d3dd1de4de6a
bd_model_1 = @reaction_network bdmodel begin
	b, 0 --> X
	d, X --> 0
end b d

# ╔═╡ ba96bd32-dc39-4a95-86c7-aeb00abf00e3
md"
It can be declared using a more concise notation:
"

# ╔═╡ eea4866b-4362-4dd7-8388-99c1bcf31b91
bd_model_2 = @reaction_network bdmodel begin
	(b,d), 0 <--> X
end b d;

# ╔═╡ f833c3a7-7528-48fb-8690-9f2c33acdaeb
md"
In the second declaration, the reaction rate is a tuple. The tuple's first element denotes the forward rate, and the second element the backward rate. The  resulting two reaction networks are identical:
"

# ╔═╡ 86785f28-8e19-4003-9a0e-f8598a187610
bd_model_1 == bd_model_2

# ╔═╡ af9591d6-5084-4fc8-97a5-1cd3eb0cfdc6
md"
## Similar reactions
Many CRNs contain a large number of reactions that differs only in the products or substrates, these can easily be bundled.
"

# ╔═╡ d20a7c6f-a4c2-444f-9f96-2e89e567f7bb
md"
Here, we have a birth-death model tracking two species (X and Y).
"

# ╔═╡ b4a04318-7e91-4b73-841b-dad3ba371cd8
xy_bd_model = @reaction_network begin
	b, 0 --> X
	b, 0 --> Y
	d, X --> 0
	d, Y --> 0
end b d

# ╔═╡ 7a6dbc82-0209-4b46-b4f4-5db7392bc34f
md"
Both the birth, and death, reactions can be bundled. Again, tuples are used to indicate bundling of reactions.
"

# ╔═╡ e2de2cad-8106-47c3-87bb-8b0c89400ebe
xy_bd_model_2 = @reaction_network begin
	b, 0 --> (X,Y)
	d, (X,Y) --> 0
end b d;

# ╔═╡ 936a82d4-e389-4f84-a757-7a1da83301fb
md"
If, e.g., the birth reactions have separate rates, this can be indicated by splitting the corresponding rate into a tuple. Thus, these two models are identical:
"

# ╔═╡ 67d58669-6183-42d4-83e8-00eda90f1c85
begin
	xy_bd_model_3 = @reaction_network begin
	bX, 0 --> X
	bY, 0 --> Y
	d, X --> 0
	d, Y --> 0
	end bX bY d

	xy_bd_model_4 = @reaction_network begin
		(bX,bY), 0 --> (X,Y)
		d, (X,Y) --> 0
	end bX bY d
end;

# ╔═╡ 7faf24a2-ff8f-4056-9cc5-a30ba267c146
md"
Bundling can be combined with reversible reactions, however, this requires nesting of tuples:
"

# ╔═╡ 505385ea-0122-4ba1-9af8-9861fba78408
xy_bd_model_5 = @reaction_network begin
	((bX,bY),d), 0 <--> (X,Y)
end bX bY d;

# ╔═╡ 973207d9-583d-4157-8e1e-2d5653053eff
md"
# Using Non-constant Reaction Rate
Reaction rates need not be a single parameter, but can be any function of parameters, species, and time.
"

# ╔═╡ b543de71-2883-468d-9a95-9eb613a775d6
md"
Here, a reaction rate can either be an expression of several parameters, or a constant.
"

# ╔═╡ b393d7ac-92d7-4b9f-9ad3-612af0357d55
bd_model_different_rates = @reaction_network begin
	b/d, 0 --> X
	1, X --> 0
end b d;

# ╔═╡ 0b3ebdbe-acdf-47a4-a3bd-a079a8d5f39e
md"
Functions can also be used:
"

# ╔═╡ 27607ab3-0a4a-4163-8e0c-0b1df2ff9503
begin
	birth(b) = b^2
	bd_model_function_rates = @reaction_network begin
		birth(b), 0 --> X
		d, X --> 0
	end b d
end;

# ╔═╡ a6b0736c-5adf-4c80-97d5-4bfae766c7c0
show_odes(bd_model_function_rates)   

# ╔═╡ b5153ab7-1e71-4ce2-8598-831c971df361
md"
## Reaction rates depending on species concentrations
E.g. where one species catalyzes a reaction, reaction rates commonly depend on species concentrations. Such reactions are possible in Catalyst.
"

# ╔═╡ e7f19045-d570-4f0d-aa87-6c663aeed353
md"
In the logistic growth model, births are not constant, but depend on species concentration. Deaths are not linear, but quadratic (indicating a carrying capacity).
"

# ╔═╡ 89d98e58-4e37-4e35-aadc-c75588213f20
begin
	logistic_growth_model = @reaction_network begin
		X, 0 --> X
		d*X, X --> 0
	end d
	show_odes(logistic_growth_model)   
end

# ╔═╡ f7cb85e6-e336-487b-bbbd-799169102b9f
md"
A biochemical reaction often requires the presence of a Catalyst, which is also a biochemical species:
"

# ╔═╡ 46239ac5-9981-4e71-9723-73e44d203dca
catalystic_conversion_model = @reaction_network begin
	k*C, S --> P
end k;

# ╔═╡ 40dca2fe-500a-4b8a-8ba1-524fafadc5cc
md"
Note that here, the Catalyst (*C*), is a species of the system, which is neither produced nor degraded:
"

# ╔═╡ b4f94075-900c-478c-98f2-5dcbe366c5da
show_odes(catalystic_conversion_model)   

# ╔═╡ 561aca56-33ac-4626-b71d-7b1790206e82
md"
Alternatively, (if *C* should be a species), it could be a part of the reaction, but with its net stoichiometry being *0*:
"

# ╔═╡ 2f74f398-520b-4b6e-ae35-189e7c973a52
catalystic_conversion_model_2 = @reaction_network begin
	k, S + C --> P + C
end k;

# ╔═╡ cdc08054-8843-44e4-8f2b-2b2b2e553d7f
show_odes(catalystic_conversion_model_2)   

# ╔═╡ 99f750e6-955f-4153-aa72-a697b0115e2a
md"
This model can also be declared with *C* as a parameter:
"

# ╔═╡ 2d8ca198-cdec-41d7-97da-391af3585d8a
catalystic_conversion_model_Cpar = @reaction_network begin
	k*C, S --> P
end k C;

# ╔═╡ 25fc399b-c082-4b89-b254-98a016406e55
show_odes(catalystic_conversion_model_Cpar)

# ╔═╡ 168d71b9-a68a-493a-8297-58cb3f1d1c78
md"
It is important to note that parameters need to be listed at the end of the system. If one forgets to list one, this will cause an error:
"

# ╔═╡ ce70f930-7e16-46ee-8e7c-628296387010
begin
	bad_model = @reaction_network begin
		(kX,kY), X <--> Y
	end kX;
	show_odes(bad_model)   
end

# ╔═╡ 29da9119-2afa-40df-8e08-0c2e9fbbc8fc
md"
(future update of Catalyst might not require this, with parameters being auto-detected)
"

# ╔═╡ bc967d1a-3250-4e26-925c-48a8fb31325b
md"
A common, non-constant, reaction rate is the Michaelis-Menten function, where the reaction rate saturates for a high concentration of the activating species. Here is an example where a single species activates its own production:
"

# ╔═╡ 48ab4be5-aaa6-4029-a547-b613d628274e
self_activation_model = @reaction_network begin
	v*X/(X+K), 0 --> X
	d, X --> 0
end v K d

# ╔═╡ 196f3459-f497-4f02-b082-b664b2c6d755
show_odes(self_activation_model)  

# ╔═╡ 4505007d-d979-4e6c-a1f3-2bf0895b09b5
md"
Since this expression is so commonly used, it is already included in Catalyst, using the mm() function.
"

# ╔═╡ 0bc1b485-6a14-47a7-8e27-92553f6373f8
self_activation_model_2 = @reaction_network begin
	mm(X,v,K), 0 --> X
	d, X --> 0
end v K d;

# ╔═╡ 90fd40b5-2126-4933-8d2f-7d9f0d1a4b97
md"
The nonlinear version of the Michaelis-Menten function, the Hill function is also available. These two models are identical:
"

# ╔═╡ 0dce9985-7a36-45b5-a7b1-3b1a29205974
begin
	self_activation_hill_model_1 = @reaction_network begin
		(v*X^n/(X^n+K^n),d), 0 <--> X
	end v K n d
	
	self_activation_hill_model_2 = @reaction_network begin
		(hill(X,v,K,n),d), 0 <--> X
	end v K n d 
end;

# ╔═╡ 16df399c-3c35-4a62-92d3-ac7f873e7d34
md"
## Reaction rates depending on time
The symbol t, when used within the Catalyst DSL, is reserved for the time variable. It can be used to make reaction rates time-dependent.
"

# ╔═╡ 8f67a2fa-6e8c-467f-8055-fcb6cde0b358
time_birth_death_model = @reaction_network begin
	b/t, 0 --> X
	d, X --> 0
end b d;

# ╔═╡ 6b6477dd-7fa2-4332-91f5-08a8d0dc8055
show_odes(time_birth_death_model)

# ╔═╡ e5231880-99ed-4d2e-939e-f78da6f1f2a2
md"
This can be used to make reactions that only occur after a certain timepoint:
"

# ╔═╡ 971da382-3cfb-4e03-8891-e8ffe821c536
timed_birth_death_model = @reaction_network begin
	b*(sign(t-t0)+1)/2, 0 --> X
	d, X --> 0
end t0 b d;

# ╔═╡ 66cfc5a9-ed0c-433d-9cde-74f031b5857f
md"
A better approach, however, is to use callbacks (which will be discussed later in the workshop). In addition, the discontinuity introduced by the step may cause problems for numeric solvers.
"

# ╔═╡ 87d6830e-53db-4aa9-be95-4275cd812294
md"
# General use of symbols in the DSL
Most Unicode characters can be used in the DSL, which can be convenient.
"

# ╔═╡ 11b3240d-74fb-40ab-a72a-c0e02ddd503a
md"
Most Unicode arrows can be used instead of \"-->\", \"<--\", and \"<-->\":
"

# ╔═╡ ec7064c2-aa8c-4496-aff6-933802e9374b
bd_model_unicode = @reaction_network begin
	(b,d), 0 ↔ X
end b d;

# ╔═╡ f412681f-eb55-4e27-a48c-8d095bff3be0
md"
This is a model of the sigmaV stress response network in *Bacillus subtilis*, where Unicode characters are extensively used to make the model declaration more comprehensible. 
"

# ╔═╡ a88b5284-1e7b-42ab-8a12-3b535af402df
σV_model = @reaction_network begin
    v₀ + hill(σᵛ,v,K,n), ∅ → (σᵛ+A)
    d, (σᵛ,A,Aσᵛ) → ∅
    (kB,kD), A + σᵛ ↔ Aσᵛ
    L*kₗ, Aσᵛ → σᵛ
end v₀ v K n kD kB kₗ d L;

# ╔═╡ 44ed2602-67d8-483a-92e4-541a3577d332
md"
# Bypassing the law of mass-action
It is likely rare that you would want to use this option, but it permits much additional freedom for declaring models.
"

# ╔═╡ 232099b8-08f9-4afe-8d9f-2e6c8cc585f0
md"
By using \"unfilled\" arrows (⇒, ⇐, ⇔), the law of mass action in an individual reaction can be bypassed. This allows reaction with arbitrary reaction rates.
"

# ╔═╡ 6a624b9b-4718-48d5-a2eb-20c0ac54f080
non_lma_model = @reaction_network begin
	p1, 0 ⇒ Y1
	p2, X ⇒ Y2
	p3, 2X ⇒ Y3
	p4, 3X ⇒ Y4
	p5, 4X + 2Z ⇒ Y5
end p1 p2 p3 p4 p5;

# ╔═╡ efa506b7-44c6-4405-92fa-140c426990cd
show_odes(non_lma_model)

# ╔═╡ dc9bbd45-b2ed-4d6a-a664-fdf0f5476517
md"
If you are using this feature extensively, it is likely that there is a better way of implementing your models.
"

# ╔═╡ 3f86541e-0307-4c3b-8bf2-65160b7a901a
md"
## Acessing the reaction system object
The \"@reaction_network\" macro creaets a *ReactionSystem* structure, which we can investigate.
"

# ╔═╡ 68797a0a-59c3-4417-9758-19ea26e97a56
binding_model = @reaction_network begin
	(p,d), 0 <--> (X,Y)
	(kB,kD), X + Y <--> XY
end p d kB kD

# ╔═╡ 83b55620-366d-45db-9d9a-e81f45cd0525
md"
Species can be fetched using the \"species\" function:
"

# ╔═╡ 65ee418e-be95-4f82-a245-e5479fa1375f
species(binding_model)

# ╔═╡ 89c3132a-c4a1-45e9-8dbe-3122258a6f43
md"
Parameters can be fetched using the \"parameters\" function:
"

# ╔═╡ e2e04178-f99e-4093-a0e0-4d3d602547d7
parameters(binding_model)

# ╔═╡ a9861c1a-6e82-4033-8d9d-31aa68290c1b
md"
Reactions can be fetched using the \"reactions\" function:
"

# ╔═╡ 7ba65dc5-1500-494e-b893-2a3ae75a642a
reactions(binding_model)

# ╔═╡ 612fe079-5d2b-4244-95f8-6e67f93a542b
md"
# Excersise: Implement networks using the Catalyst DSL
Implement the following reaction networks using the Catalyst DSL. You can reveal the network box to see the network. If you use the same name for your model (optional input between \"@reaction_network\" and \"begin\"), you can check equality using \"==\". Try writing the network by combining similar reactions.
"

# ╔═╡ 66728203-e607-4b03-8520-85c2864b89ec
excersise_rn_1 = @reaction_network ern1 begin
	p, 0 --> (X,Y)
	d, (X,Y) --> 0
	(kB,kD), X + Y <--> XY
end p d kB kD

# ╔═╡ e75d0780-28f9-4578-8cad-f52736cd2160
# answer_rn_1 = @reaction_network ern1 begin
# 	 # Write network here.
# end 
# answer_rn_1 == excersise_rn_1

# ╔═╡ 65f61e36-7794-41c2-b04f-1422ffd45870
excersise_rn_2 = @reaction_network ern2 begin
	b, S + I --> 2I
	k, I --> R
end b k

# ╔═╡ 0a0144cb-25bd-413b-adbd-31436b3c7cd4


# ╔═╡ 8899e7dd-8635-45d6-9e96-c86d6a0c91b0
excersise_rn_3 = @reaction_network ern3 begin
	(p,d), 0 <--> C
	mm(C,v,K), S --> P
end p d v K

# ╔═╡ bdb26875-46c2-4398-a3dc-1a462d7dd889


# ╔═╡ 6c4881a5-f6ba-45cb-aaa4-376bb78e8050
excersise_rn_4 = @reaction_network ern4 begin
	hill(K,v,Z,n), 0 --> X
	hill(K,v,X,n), 0 --> Y
	hill(K,v,Y,n), 0 --> Z
	d, 0 --> (X,Y,Z)
end v K n d

# ╔═╡ bf7922f9-a721-4dc0-a766-3bc005927936


# ╔═╡ Cell order:
# ╟─37a68fe9-328d-4e10-84da-113dae3b1a9e
# ╟─3283a80a-09cc-11ed-2959-6525dd49b3fd
# ╟─6ec927cc-a67b-44ba-b750-6947a1518dbc
# ╟─76af30dd-f995-4a48-b922-76c2e84b1693
# ╟─073c9ca2-7504-42ec-ad9f-e9e3014303f3
# ╟─35be8093-beb7-45f9-bfa1-51e91e27292a
# ╟─0879683c-b797-4ba1-8e86-8a16006bbdc5
# ╟─c78dbeaf-bda0-4be7-816b-553bf7300829
# ╟─344dbafb-0d22-40f4-a250-e16ac366a41a
# ╟─66e4e1c9-3678-4235-8bc6-71c24b1cd1d5
# ╠═32ec973b-4e14-461a-bf30-7262265c0fe0
# ╟─06d766b4-261b-46ac-a948-ea51ffe42572
# ╠═a03aec3b-82c6-4dc6-b046-20f7f7b04c5f
# ╟─fb29bcf8-e30c-4ff5-bed5-82314d371204
# ╟─620f072c-1d4b-4ae3-895e-c792471502fe
# ╠═c0ab7ba6-0437-44cc-85d9-af6d2753339a
# ╟─83b8f3fd-bff1-4588-b951-87863ff885b5
# ╟─024163ae-0ad4-41fd-a48b-0ed85876d960
# ╟─2feb2fe7-2c07-4d7e-9b19-e8e7d68f1e24
# ╠═b838469c-7658-41a5-af1c-272fe1ceff68
# ╟─fbdcc69a-6849-48c4-8e8c-de3c31bb77c2
# ╟─3f1c0b3c-f988-4180-b637-ba6e11f8ada8
# ╟─7cfadb90-4283-4471-b707-373b01aa4fd8
# ╟─ad17c467-0852-4414-a2c4-df8b36329846
# ╠═a4b6ad01-2713-4927-89c2-d3dd1de4de6a
# ╟─ba96bd32-dc39-4a95-86c7-aeb00abf00e3
# ╠═eea4866b-4362-4dd7-8388-99c1bcf31b91
# ╟─f833c3a7-7528-48fb-8690-9f2c33acdaeb
# ╠═86785f28-8e19-4003-9a0e-f8598a187610
# ╟─af9591d6-5084-4fc8-97a5-1cd3eb0cfdc6
# ╟─d20a7c6f-a4c2-444f-9f96-2e89e567f7bb
# ╠═b4a04318-7e91-4b73-841b-dad3ba371cd8
# ╟─7a6dbc82-0209-4b46-b4f4-5db7392bc34f
# ╠═e2de2cad-8106-47c3-87bb-8b0c89400ebe
# ╟─936a82d4-e389-4f84-a757-7a1da83301fb
# ╠═67d58669-6183-42d4-83e8-00eda90f1c85
# ╟─7faf24a2-ff8f-4056-9cc5-a30ba267c146
# ╠═505385ea-0122-4ba1-9af8-9861fba78408
# ╟─973207d9-583d-4157-8e1e-2d5653053eff
# ╟─b543de71-2883-468d-9a95-9eb613a775d6
# ╠═b393d7ac-92d7-4b9f-9ad3-612af0357d55
# ╟─0b3ebdbe-acdf-47a4-a3bd-a079a8d5f39e
# ╠═27607ab3-0a4a-4163-8e0c-0b1df2ff9503
# ╟─a6b0736c-5adf-4c80-97d5-4bfae766c7c0
# ╟─b5153ab7-1e71-4ce2-8598-831c971df361
# ╟─e7f19045-d570-4f0d-aa87-6c663aeed353
# ╠═89d98e58-4e37-4e35-aadc-c75588213f20
# ╟─f7cb85e6-e336-487b-bbbd-799169102b9f
# ╠═46239ac5-9981-4e71-9723-73e44d203dca
# ╟─40dca2fe-500a-4b8a-8ba1-524fafadc5cc
# ╟─b4f94075-900c-478c-98f2-5dcbe366c5da
# ╟─561aca56-33ac-4626-b71d-7b1790206e82
# ╠═2f74f398-520b-4b6e-ae35-189e7c973a52
# ╟─cdc08054-8843-44e4-8f2b-2b2b2e553d7f
# ╟─99f750e6-955f-4153-aa72-a697b0115e2a
# ╠═2d8ca198-cdec-41d7-97da-391af3585d8a
# ╟─25fc399b-c082-4b89-b254-98a016406e55
# ╟─168d71b9-a68a-493a-8297-58cb3f1d1c78
# ╠═ce70f930-7e16-46ee-8e7c-628296387010
# ╟─29da9119-2afa-40df-8e08-0c2e9fbbc8fc
# ╟─bc967d1a-3250-4e26-925c-48a8fb31325b
# ╠═48ab4be5-aaa6-4029-a547-b613d628274e
# ╟─196f3459-f497-4f02-b082-b664b2c6d755
# ╟─4505007d-d979-4e6c-a1f3-2bf0895b09b5
# ╠═0bc1b485-6a14-47a7-8e27-92553f6373f8
# ╟─90fd40b5-2126-4933-8d2f-7d9f0d1a4b97
# ╠═0dce9985-7a36-45b5-a7b1-3b1a29205974
# ╟─16df399c-3c35-4a62-92d3-ac7f873e7d34
# ╠═8f67a2fa-6e8c-467f-8055-fcb6cde0b358
# ╟─6b6477dd-7fa2-4332-91f5-08a8d0dc8055
# ╟─e5231880-99ed-4d2e-939e-f78da6f1f2a2
# ╠═971da382-3cfb-4e03-8891-e8ffe821c536
# ╟─66cfc5a9-ed0c-433d-9cde-74f031b5857f
# ╟─87d6830e-53db-4aa9-be95-4275cd812294
# ╟─11b3240d-74fb-40ab-a72a-c0e02ddd503a
# ╠═ec7064c2-aa8c-4496-aff6-933802e9374b
# ╟─f412681f-eb55-4e27-a48c-8d095bff3be0
# ╠═a88b5284-1e7b-42ab-8a12-3b535af402df
# ╟─44ed2602-67d8-483a-92e4-541a3577d332
# ╟─232099b8-08f9-4afe-8d9f-2e6c8cc585f0
# ╠═6a624b9b-4718-48d5-a2eb-20c0ac54f080
# ╠═efa506b7-44c6-4405-92fa-140c426990cd
# ╟─dc9bbd45-b2ed-4d6a-a664-fdf0f5476517
# ╟─3f86541e-0307-4c3b-8bf2-65160b7a901a
# ╠═68797a0a-59c3-4417-9758-19ea26e97a56
# ╠═83b55620-366d-45db-9d9a-e81f45cd0525
# ╠═65ee418e-be95-4f82-a245-e5479fa1375f
# ╟─89c3132a-c4a1-45e9-8dbe-3122258a6f43
# ╠═e2e04178-f99e-4093-a0e0-4d3d602547d7
# ╟─a9861c1a-6e82-4033-8d9d-31aa68290c1b
# ╠═7ba65dc5-1500-494e-b893-2a3ae75a642a
# ╟─612fe079-5d2b-4244-95f8-6e67f93a542b
# ╟─66728203-e607-4b03-8520-85c2864b89ec
# ╠═e75d0780-28f9-4578-8cad-f52736cd2160
# ╟─65f61e36-7794-41c2-b04f-1422ffd45870
# ╠═0a0144cb-25bd-413b-adbd-31436b3c7cd4
# ╟─8899e7dd-8635-45d6-9e96-c86d6a0c91b0
# ╠═bdb26875-46c2-4398-a3dc-1a462d7dd889
# ╟─6c4881a5-f6ba-45cb-aaa4-376bb78e8050
# ╠═bf7922f9-a721-4dc0-a766-3bc005927936
