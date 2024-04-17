# tweedieCARpanel

This is a bare-bones (really, just one function) package that I developed for a specific use-case during a consulting project. I thought the source code might be useful to others, so it is now public.

The function fits a Tweedie GLM with a latent spatial conditional autoregressive model based on a spatial adjacency matrix. The spatial diffusion effect can be permanent or vary in time by changing the `fixed_speff` argument to `FALSE` and specifying the number of time steps.

**Note**: AS OF 17 APRIL 2024, I STILL NEED TO ADD CODE FOR THE TIME-VARYING SPATIAL DIFFUSION MODEL!

# Installation

To download the package with precompiled C++ code, use:

`devtools::install_github("Dusty-Gannon/tweedieCARpanel", build = FALSE)`

