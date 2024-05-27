# GalerkinAppr
Galerkin approximator for 2D elasticity problems in solid mechanics

This program is intended for educational purposes. In Galerkin approximations the approximation functions for displacements are defined in the whole domain. This is not the case in popular finite element methods (which are in fact based on Galerkin's approximation), the main difference being the use of element-wise defined functions.

Bare Galerkin approximations are feasible for simple geometries only, but the process of choosing adequate functions for the problem at hand has a high educational value. This program is intended to make easier the process of trial-and-test the proposed functions. Usually, we want to use no more than a handful of functions but, even in this case, operating the approximation by hand would be a discouraging scenario.

The program runs on Python3, and requires Tkinter, MatPlotLib, NumPy and SymPy libraries.
To get started please run the pre-built example from the interface. The explanations in the popup window are the closest to a 'user manual' available.

I hope you'll find GalerkinAppr interesting.
The author.
JC del Caño
_____________________
