# PressureFieldContact.jl

[![Build Status](https://travis-ci.org/ryanelandt/PressureFieldContact.jl.svg?branch=master)](https://travis-ci.org/ryanelandt/PressureFieldContact.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/PressureFieldContact.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/PressureFieldContact.jl?branch=master)

This module implements the elastic foundation contact model for rigid body dynamics described in [this video](https://drive.google.com/open?id=1R_q9eIaIBnTLhvTE5U2uzUsbZOM8hdeV).
[This paper](https://arxiv.org/pdf/1904.11433.pdf) describes the method in greater detail.
The surface of bodies are represented with a triangular mesh.
The compliant portion of bodies is represented with a tetrahedral mesh.
The easiest way to get started is to run the boxes.jl example in the test directory.
More examples will come soon.

### Friction Models

This package can model friction using either a regularized Coulomb friction model or a bristle friction model.
A regularized Coulomb friction model is simple at the cost of allowing creep (i.e. a block will slide down a ramp irrespective of slope).
The bristle friction model adds extra state variables that model deformation.
This friction model allows objects to have zero tangential velocity even when applied tangential forces are not zero.
Try both to see which is best for your application.
