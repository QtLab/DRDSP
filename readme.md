DRDSP
=====

Dimensionality Reduction for Dynamical Systems with Parameters

Introduction
------------

DRDSP is a C++ library for producing a low-dimensional dynamical system that models the long-term behaviour (attractors) of a given high-dimensional system. Systems with parameter spaces -- resulting in families of attractors -- can also be reduced.

This library implements the approach described in [this paper](http://dx.doi.org/10.1137/130913675).

The steps involved in the method are:

1. Generate some data from the attractors.
2. Compute the secants from the data sets.
3. Use the secants to find a projection.
4. Use the projection, the data, and the original model to compute some reduced data.
5. Use the reduced data to obtain the reduced model.

The method is compatible with any autonomous smooth deterministic dynamical system corresponding to the flow of a vector field.

Non-Euclidean state spaces are supported (e.g. systems with angular state variables).

Non-autonomous systems with periodic time-dependence are supported by encapsulating the explicit time-dependence as additional state (resulting in an autonomous system).

Usage
-----

See [the wiki](https://github.com/cwzx/DRDSP/wiki).

For more detailed reference see `docs/`.


Examples
--------

Examples are provided in `examples/`.

* `examples/rossler/` -- The Rossler system. 3-dimensional state space.
* `examples/pendulum/` -- A double pendulum with friction and external forcing. 5-dimensional non-Euclidean state space.
* `examples/brusselator/` -- The Brusselator PDE on a 2D physical space.
* `examples/kuramoto/` -- The Kuramoto model. Coupled oscillators with external forcing.
* `examples/goodfellow/` -- A network of nonlinear oscillators.
* `examples/dynamo/` -- A very complicated discretized PDE, with a 13120-dimensional state space.


Extra Features
--------------

* `HistogramGenerator` -- Generates a histogram from a data set.

* `BifurcationDiagramGenerator` -- Performs simulations for multiple parameter values and samples values to be plotted on a bifurcation diagram. The sampling condition and value functions must be specified. A custom integrator may also be specified (RK4 by default).


Build
-----

The `include` directory contains the header files that should be copied to your include path.

The `build` directory contains project files for Visual Studio 2013, which will build a static library to `lib/` and executables for the examples to `bin/`. The static library should be linked in your build process.

At the moment there are no build files for other compilers/platforms. You will find the source for the library in `src/` if you wish to make your own.

### Dependencies and Requirements

DRDSP uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), a C++ template library for linear algebra. This is included in the repository.

A number of C++11 features are used by the library, and maybe one or two C++14 features.


Acknowledgements
----------------

The development of this library was supported by EPSRC grant EP/G026238/1 as part of the myGrid project.

