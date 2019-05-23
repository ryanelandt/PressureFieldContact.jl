# Geometry

This contact implementation represents rigid objects with triangular meshes and compliant objects with tetrahedral meshes.
This package uses the custom geometry data type `eMesh`.
Convenience functions to create meshes for common shapes are documented below.

## `eMesh` Data Type

```@autodocs
Modules = [PressureFieldContact.Geometry]
Order   = [:type]
Pages   = ["mesh.jl"]
```

## Simple Shapes

```@docs
eMesh_box
eMesh_sphere
eMesh_half_plane
eMesh_cylinder
```

## Swept Meshes

```@docs
create_swept_mesh
f_swept_triv
f_swept_circle
```

## Rotationally Symmetric Meshes

```@docs
```

## Mesh from an STL file

For an example of this this example:


## Transforming Meshes

```@docs
transform!
```

```@docs
basic_dh
```
