DRDSP
=====

Dimensionality Reduction for Dynamical Systems with Parameters

Intro
-----

DRDSP is a C++ library for producing a low-dimensional dynamical system that models the long-term behaviour of a given high-dimensional system.

The steps involved in the method are:

* Generate some data from the attractors.
* Compute the secants from the data sets.
* Use the data and its secants to find a projection.
* Use the projection and the data to reduce the data.
* Use the reduced data and the original model to obtain the reduced model.

Usage
-----

### Interfaces to be implemented

In order to use the library with your example system, you need to implement some interfaces to specify the dynamics of the system.

* `DRDSP/dynamics/embedding.h`
** `Embedding` -- An embedding of the state space into R^n. The default implementation is the identity map, which is for systems whose state space is already R^n. For example, if your model contains angular variables, these must be embedded.

* `DRDSP/dynamics/dynamicalSystem.h`
** `WrapFunction` -- If your system contains angular variables, they must be wrapped to the range [-pi,pi] to avoid a loss of precision. The `WrapFunction` class is a function object that takes a state and performs the wrapping. The default implementation leaves the state unchanged.

* `DRDSP/dynamics/model.h`
** `Model*` -- This is where the dynamics of the original model are specified. Both the vector field and its partial derivatives (Jacobian) must be implemented.

*** `Model` -- A single vector field, no parameters.
*** `ModelCW` -- A single vector field, no parameters. Evaluate component-wise.
*** `ModelParameterized` -- A family of vector fields with parameters.
*** `ModelParameterizedCW` -- A family of vector fields with parameters. Evaluate component-wise.

### Data sets

The method requires data samples from the attractors of the original model. These data sets can either be loaded in from a file, or generated at run-time.

To generate a data set, the `DataGenerator' class will perform a simulation of the original model using an RK4 integrator for multiple parameter values.

Once data has been obtained, the secants are computed using the `Secants` class.

Secants can be culled using one of the `Secants::CullSecants*` methods.

### Secant-based Projection

The `ProjSecant` class determines a projection from a set of secants. The projection is encoded in an orthonormal matrix, `ProjSecant::W`.

### Computing Reduced Data
`ReducedData` and `ReducedDataSystem` classes take the original model, the data sets, and an orthonormal matrix and compute the necessary data required to find the reduced model.

### Find Reduced Model

The `ModelRBFProduced` class takes a `ReducedData` and produces a `ModelRBF`.
The `ModelReducedProducer` class takes a `ReducedDataSystem` and produces a `ModelReduced`.


Examples
--------

Examples are provided in `examples/`.

* `examples/pendulum/` -- A double pendulum with friction and external forcing.
* `examples/brusselator/` -- The Brusselator PDE on a 2D physical space with a 64x64 discretization.
* `examples/ks/` -- The Kuramoto model. Coupled oscillators with external forcing.

Build
-----

The `include` directory contains the header files that should be copied to your include path.

The `build` directory contains project files for Visual Studio 2012, which will build a static library to `lib/` and executables for the examples to `bin/`. The static library should be linked in your build process.


