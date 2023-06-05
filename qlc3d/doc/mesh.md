# Mesh and Geometry


## Geometry
### Material Names
Every region of the modelling window (volume, surface) must have a "material" associated
with it. This is used to identify e.g. LC regions, Dielectric regions, boundary conditions etc.

The supported materials are defined in the file `material_number.h` file and are currently:

#### Volume materials:

- `Domain1` - LC region 
- `Dielectric1` - `Dielectric7` - Dielectric regions. Each can have a different dielectric constant assigned to it. 

#### Surface materials:
- `Periodic` - Periodic boundary condition. Note that the actual mesh must also be periodic in these regions.
- `Neuman` Neumann boundary condition for potential calculation and "free surface" for LC director calculation.
- `Electrode1` - `Electrode9` - Electrode regions. Each can have a different voltage assigned to it.
- `FixLC1` - `FixLC9` - LC director anchoring surfaces. Each FixLC surface can have different anchoring parameters assigned to it. 

Note that it is possible to combine Electrode and FixLC materials such that one surface has both a voltage and a fixed LC director.
How this is done depends on the meshing software used to create the mesh. See the included examples for more details.

### Periodicity
In general, a modelling window of any shape is supported, however axis-aligned cuboid windows are the most common in LC device modelling.
When using am axis-aligned cuboid mesh, periodic boundary conditions are supported in the following directions:

#### Periodic ymin-ymax boundaries
When the modelling window is periodic along a single axis, the axis must be y. This can for example be 
used to model a 2D slice in the x-z plane through a 3D device, where the LC orientation is effectively uniform along the y-axis. 
The xmin and xmax boundaries can then be set to e.g. Neumann boundary conditions while the zmin and zmax boudaries can 
represent LC alignment layers or electrodes.

#### Periodic ymin-ymax, xmin-xmax boundaries

When two axes of periodicity exist, they must be along the x and y axes.

This can be an extension of the 2D slice example above, but with the Neumann boundaries replaced by periodic boundaries.
Alternatively this can be used to represent a periodically repeating section of a device, e.g. one or more pixels in a display. 

#### Periodic ymin-ymax, xmin-xmax, zmin-zmax boundaries.
When all three axes are periodic, the modelling window typically approximates a smaller region of a "bulk" LC material which 
isn't in contact with any electrodes or alignment layers. It may however contain volume regions of dielectric material, such as
"particles" floating withing the LC material. 

## Mesh
- First order, linear, tetrahedral and triangular elements are used
- support for Gmsh ASCII mesh format 4.1, generated using the free Gmsh software (http://gmsh.info/). For example, 
Gmsh version 4.8.4, included in Ubuntu 22.04 has been verified to work. 

