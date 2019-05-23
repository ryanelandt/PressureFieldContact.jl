# Add Contact/Friction

This page explains how to take geometry and as `eMesh` and make it engage in contact.
This process has two steps: 1.) add the `eMesh` to the contact geometries and 2.) add friction pairs to contact geometry pairs.

## Add Contact

```@docs
add_body_contact!
add_contact!
add_body!
as_tri_eMesh
as_tet_eMesh
```

## Add Friction

Friction applies between user-selected contact pairs.
This model has two friction contact models.
The first is a piece-wise linear regularization of Coulomb friction.
The second is a 6-DOF bristle friction model.

```@docs
add_friction_regularize!
add_friction_bristle!
```
