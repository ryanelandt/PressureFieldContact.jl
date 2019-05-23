# PressureFieldContact.jl

[![Build Status](https://travis-ci.org/ryanelandt/PressureFieldContact.jl.svg?branch=master)](https://travis-ci.org/ryanelandt/PressureFieldContact.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/PressureFieldContact.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/PressureFieldContact.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://ryanelandt.github.io/PressureFieldContact.jl/dev)

This module implements the elastic foundation-themed contact model for rigid body dynamics described in [this short video](https://drive.google.com/open?id=1R_q9eIaIBnTLhvTE5U2uzUsbZOM8hdeV).
[This paper](https://arxiv.org/pdf/1904.11433.pdf) describes the method in greater detail.
UNDER CONSTRUCTION: See the latest [documentation](https://ryanelandt.github.io/PressureFieldContact.jl/dev) for installation instructions, a quick-start guide and summary of how the different pieces of this method work.

### Summary

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
