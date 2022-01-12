# heom-lab
A Matlab implementation of the Hierarchical Equations of Motion (HEOM) method for modelling open system quantum dynamics.

## Summary of heom-lab
heom-lab provides a matlab implementation of the HEOM method including some non-standard features. Such features include:
* Dynamics for arbitrary numbers of baths with arbitrary system coupling operators.
* Dynamics with Debye spectral densities as well as underdamped and overdamped Brownian oscillator baths.
* Integration of the equations of motion with an adaptive short iterative Arnoldi (SIA) integrator.
* Truncation of the hierarchy with the standard level-based truncation, as well as more efficient truncation schemes for larger reorganisation energy/high temperature problems based on frequency cut-off of the hierarchy and a coupling/frequency weighted scheme.
