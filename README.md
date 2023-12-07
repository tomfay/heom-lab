# heom-lab
A Matlab implementation of the Hierarchical Equations of Motion (HEOM) method for modelling open system quantum dynamics.

_A note of caution_: This code is very much a work in progress. Any reports of bugs and issues, or requests for features are welcome. If you're interested in using the code or helping to develop it, please get in touch at tom.patrick.fay AT gmail.com.

## Summary of heom-lab
heom-lab provides a matlab implementation of the HEOM method including some non-standard features. Such features include:
* Dynamics for arbitrary numbers of baths with arbitrary system coupling operators.
* Dynamics with Debye spectral densities as well as underdamped and overdamped Brownian oscillator baths.
* Integration of the equations of motion with an adaptive short iterative Arnoldi (SIA) integrator.
* Truncation of the hierarchy with the standard level-based truncation, as well as more efficient truncation schemes for larger reorganisation energy/high temperature problems based on frequency cut-off of the hierarchy and a coupling/frequency weighted scheme.

## Citations
If you do use this code for results in any publications, I'd appreciate it if you could cite: 
"A simple improved low temperature correction for the hierarchical equations of motion" Thomas P. Fay, J. Chem. Phys. 157, 054108 (2022) https://doi.org/10.1063/5.0100365
If you use the hybrid-HEOM methods, please also cite:
"Coupled charge and energy transfer dynamics in light harvesting complexes from a hybrid hierarchical equations of motion approach" Thomas. P Fay and David T. Limmer J. Chem. Phys. 157, 174104 (2022) https://doi.org/10.1063/5.0117659
Thanks!


