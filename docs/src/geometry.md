# Geometry

This contact implementation represents rigid objects with triangular meshes and compliant objects with tetrahedral meshes.
This package uses the custom geometry data type `eMesh`.
Convenience functions to create meshes for common shapes are documented below.

```@autodocs
Modules = [PressureFieldContact.Geometry]
Order   = [:type]
Pages   = ["mesh.jl"]
```

## Basic Shapes

```@docs
eMesh_box
eMesh_sphere
eMesh_half_plane
eMesh_cylinder
```

## Transforming Meshes

```@docs
transform!
```

## Swept Meshes



## Rotationally Symmetric Meshes



## Mesh from an STL file

For an example of this this example:
