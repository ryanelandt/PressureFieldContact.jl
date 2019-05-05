# PressureFieldContact.jl

[![Build Status](https://travis-ci.com/ryanelandt/PressureFieldContact.jl.svg?branch=master)](https://travis-ci.com/ryanelandt/PressureFieldContact.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/PressureFieldContact.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/PressureFieldContact.jl?branch=master)

This module implements an elastic foundation contact model for rigid body dynamics.
The surface of bodies is represented with a triangular mesh.
The compliant portion of bodies is represented with a tetrahedral mesh.
The easiest way to get started is to run one of the examples.

### Friction Models

This package can model friction using either a regularized Coulomb friction model or a bristle friction model.
A regularized Coulomb friction model is simple at the cost of allowing creep (i.e. a block will slide down a ramp irrespective of slope).
The bristle friction model adds extra state variables that model deformation.
This friction model allows objects to have zero tangential velocity even when applied tangential forces are not zero.
Use this model for gripping tasks.
