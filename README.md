# interFoamMod
Modified interFoam solver for handling highly viscous fluids based on the modifications suggested in https://doi.org/10.1115/1.4053548.

Difference from default interFoam code is made only with respect to introducing a smoothing method for calculation of the curvature under interfaceProperties.H, and introducing an exponential averaging equation for property estimation based on the viscosity ratios. Refer to the paper for the complete details.
