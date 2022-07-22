# JuliaCon 2022 Catalyst Workshop
This repo contains the notebooks used in the workshop. They can be installed and used to follow along as follows:

1. Click the green code button on Github and download the package to your computer. If needed, unzip the file.
2. Start Julia within the JuliaCon2022_Catalyst_Workshop folder, or navigate to this folder within Julia.
3. Enter
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Pluto
Pluto.run()
```
4. At this point Pluto should open in your web browser, and you can load any of the workshop notebooks.

Note, these tutorials use Catalyst 12.1.2.
