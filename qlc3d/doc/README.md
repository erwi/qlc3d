# Running qlc3d

It is recommended that qlc3d is run from the command line. The reason is because when qlc3d encounters problems, error messages that may help resolving the problem are printed on the screen. If qlc3d is started by simply double clicking on the executable, the opened window will close before there is time to read the error message.

In general, qlc3d is run by typing to the command prompt:

```
qlc3d <settings file name>
```

where `<settings file name>` is the name of a text file that contains the parameter values required to run the simulation. This file is called the settings file in this document, but the actual file name can be anything (e.g. “settings.txt” will do). The contents of this file are described later in this document. 

In addition to the settings file, qlc3d requires a mesh file that contains the finite elements mesh for the simulated structure. Typically this file would have been made using some external mesh generator program, such as [GiD](https://www.gidhome.com/) or [Gmsh](https://gmsh.info/). The name of the mesh file, among other things, is specified in the settings file.

## Suggested Directory Structure

It is convenient to create a single descriptively named *project directory* that contains a single settings file and a single mesh file. It is of course possible to keep multiple settings/mesh files in this directory (e.g. multiple related but somehow different meshes), but in general it is easier to keep unrelated simulations in separate directories. This ensures that all the files needed to run a simulation (or series of related simulations) and the result files output by qlc3d are kept well organised and together. 

When qlc3d is run from the project directory, by default a result directory named “res” is created within the project directory and all the result files are written there. Multiple result directories with different names can be created by specifying this in the settings file (see the `SaveDir` setting). This makes it possible to run qlc3d with different parameters without overwriting previous results. 

In addition to result files, qlc3d also makes copies of both the current settings and the current mesh files and places them in the result directory. This means that each result directory also contains all the input files used to generate the result files. This makes it easy to repeat a simulation or check the value of some parameter etc. at a later time (for example when writing up a final report, some months after the actual results were obtained). 

# Settings Files

qlc3d is controlled by a single text file that contains all the required settings and parameters. This file is called the *settings file*. 

Some or the parameters described below are mandatory (such as the name of the mesh file to use), but not all. When a parameter value is not defined in the settings file, “reasonable” default value is assumed for the parameter. When a parameter is defined, its value overrides the default. 

**The exception...**

In some cases, the values of some parameters result in other (conflicting or nonsensical) parameters being ignored, even when defined. For example, when the type of an anchoring surface is defined to be “strong” (implying infinite anchoring energy strength),  its (finite) strength value as defined in the settings file is ignored and does not affect the simulation.

## Settings File Format

A settings file is a simple text file that can be edited using any text editor (e.g. Notepad++ on widows or gedit/Kate/nano/vi(m) on Linux). The file contains a number of '='-sign separated key-value pair assignments:
```
	<key> = <value>
```

where `<key>` is the name of a parameter and `<value>` defines its value. Capitalisation of both keys and values is ignored, so that the following are all equal:
```
	key = value
	Key = Value
	KeY = VALUE
```

The values can be scalars, vectors/lists or text, depending on their meaning. Vector/list values start with ‘[‘and end in ‘]’, and each value is separated by a comma. In addition, anything on a line following the '#’ –character is interpreted as a comment and does not affect the simulation. 

Examples of valid vector, scalar and text key/value pairs with comments:

```
	# This is a comment
	SomeVector = [1, 2, 3] # This is also a comment
	SomeScalar = 1e-3
	SomeText = text
	SomeListOfText = [text1, text2, text3]
```    

The settings file is not interpreted, so writing expressions including arithmetic operations results in an error and referring to other previously defined key name as a variable is treated as a key/value pair with a text value, so it is not possible to do  scripting with variables in the settings file.

```	
	key1 	= 2	# This is a valid assignment
	key2 = 1+2	# error, arithmetic not supported
	key3 = key1 # not error, but assigns text to key3
```

### String Substitution

As explained above, qlc3d does not support scripting or arithmetic expressions in the settings file. However, as a "poor man's substitute", string-substitution from environment variables is supported. This makes it relatively easy to set up more complicated simulations that may, for example, perform sweeps of ranges of values for some settings. 

When qlc3d reads in a settings file, it first checks the file contents for specially formatted sub-strings, starting with the `$` character followed by an environment variable name enclosed in curly braces `{` and `}`, (e.g. `${SOME_ENVIRONMENT_VARIABLE_NAME}`), it then substitutes it with the value corresponding to the named environment variable.


For example, if the environment variable `MY_ELECTRODE_POTENTIAL` has been defined as `1.5`, then the line
```
    E1.Pot = [${MY_ELECTRODE_POTENTIAL}] # switching potential for electrode 1
```
will be substituted as 
```
    E1.Pot = [1.5] # switching potential for electrode 1
```
Examples of using this can be found in the included `examples` directory, where python scripts are used to repeatedly call qlc3d with different values defined at each run for some parameters.

# Index of Settings File Keys

Below, the possible key-value pairs are listed with brief explanations. These are roughly grouped by their use.

---
## General Settings
### **MeshName**
Mandatory

Defines the file name of the finite elements mesh sed. **Note**: the path is w.r.t. the current directory where qlc3d is run from.

Example 1: the mesh file "msh.msh" is in the current directory
```
	MeshName = mesh.msh
```
Example 2: the mesh file "mesh.msh" is in a sub-directory named "meshes"
```
	MeshName = meshes/mesh.msh
```

### **EndCriterion**
Optional (default = Iteration, TODO: check this)

Specifies what criterion is used to end a simulation. This parameter is used in conjunction with the `EndValue` parameter. Supported values for this parameter are:

|Value 	| Explanation	|
|---	| ---			|
|Iterations| The simulation finishes when number of iterations or steps has reached the value specified in `EndValue`|
|Time| The simulation finishes when the simulated time is greater or equal to the value specified in `EndValue`|
|Change|The simulation finishes when the Q-tensor appears to have reached a steady state. The largest change in any of the Q-tensor components between two successive iterations or steps is monitored and the value specified for EndValue is used as the threshold for determining when steady state is reached.|


When this option is used, both the current Q-tensor change and the target value are displayed on the screen at the end of each step.

**Note** that if the EndValue used is too small, it is possible that this criterion is never met. This can be caused e.g. due to a too coarse a mesh or something else preventing convergence. 


### **EndValue** ###

Optional (default value EndValue = 5 TODO:check this).

Defines the numerical threshold value used in conjunction with the `EndCriterion` parameter to determine when to end a simulation.  The meaning of this value depends on the value used for `EndCriterion`:

Example 1: Simulate for 1000 iterations or time steps:

```
	EndCriterion = Iterations
	EndValue = 1000
```

Example 2. Simulate a 15 ms time period:

```
	EndCriterion = Time
	EndValue = 15e-3
```

Example3. Simulate until largest change in Q-tensor is below 6 significant digits:

```
	EndCriterion = Change
	EndValue = 1e-6
```

### **dt** ###
Optional (default value: dt = 1e-9)

Numerical value determining the initial time step duration in units of seconds. 

qlc3d uses an adaptive time stepping scheme that automatically changes the size of the time step at each step in an attempt to maintain a constant change in the Q-tensor for each time step. That is, the size of the time step is automatically decreased when the Q-tensor varies rapidly to maintain accuracy. Similarly, the size of the time step is increased when the Q-tensor changes slowly in order to speed up the simulation. The way this algorithm works is controlled by the parameters `dtLimits`, `TargetdQ` and `dtFunction`.

If dt is set to zero, a steady-state solver mode (Newton-Raphson method) is used instead. This mode can be used to more quickly find a stable steady state solution, when the dynamics of the LC material are not of interest. The steady state mode may fail to converge to a solution for a number of reasons, including:

1. Initial starting condition is not good enough.
2. The mesh is not fine enough for the problem (e.g. in the presence of defects).

Hint:
It is often a good idea to first run a simulation for a short while (e.g. for one milliseconds) with dt > 0.  Then using a the last result file as an initial starting LC configuration and setting dt = 0 to try to find the final steady state solution.  

See also the `LoadQ`, `EndCriterion` and `EndValue` parameters to do this. 


### **dtLimits = [1e-09, 1e-3]** ###
Vector of length 2, specifying the minimum and maximum values that `dt` can take. Adaptive time stepping will adjust the value of dt during a simulation if initial `dt` is not 0.

### **TargetdQ = 1e-3** ###
When time-stepping `(dt > 0)`, the value of `dt`, is adapted after each step so that the maximum change in the Q-tensor approaches this value. I.e., if `dQ > TargetdQ`, `dt` is decreased and vice versa.

When trying to directly solve for a steady state solution using the Newton method (dt = 0), this parameter controls the degree of damping applied. The change in Q-tensor is scaled so that it is always less or equal to TargetdQ.  This helps preventing the solution from diverging. 

### **dtFunction = [0.25, 0.8, 1.2, 10]** ##
Array of length 4. The value of the current time step size `dt` is changed according to `dt(n+1) = S(R) * dt(n)`, where `R` is the ratio `R = dQ/TargetdQ`, `n` is time-step number and `S` is a piecewise linear function whose extreme values are limited between 0 and 2.

```	
	dtFunction[0] = value of R where S = 0;
	dtFunction[1] = minimum value of R where S = R;
	dtFunction[2] = maximum value of R where S = R;
	dtFunction[3] = value of R where S = 2;
```

Additionally, the limits set in `dtLimits` are enforced.

### **MaxError = 1.0000e-03** ###
numerical value used as an accuracy parameter for the nonlinear Crank-Nicholson time-stepping. Newton iterations are performed within the time step until dQ is less or equal to `MaxError`,

### **numAssemblyThreads = 0** ###
An optional numerical value, larger or equal to 0, that specifies the number of threads to use in assembling the matrix problems.  The default value is 0, which results in all the available threads being used.

### **NumMatrixSolverThreads = 0** ###
An optional numerical value, larger or equal to 0, that specifies the number of threads to use in solving the matrix problems.  The default value is 0, which results in all the available threads being used. Note that, optimum performance is achieved using a value larger than 0 but smaller than the number of cores/hardware threads available on the system.

---
## Result Input/Output ##
Settings related to specifying result file output formats are given below.


### **SaveDir = name of results directory** ###
String value determining a sub-directory where results are saved. If the directory does not exist, it will be created. If result files from previous simulations already exist in this directory, these will be overwritten. This value is empty by default, in which case results are written in the sub-directory `./res` (which is created automatically if it does not exist).


### **SaveFormat = [type1, type2 … typeX]** ###
Vector of strings specifying file format of saved result files. Possible file types are:
`LCview`, `LCviewTXT`, `RegularVTK`, `RegularVecMat`, `DirStackZ`.  When `SaveFormat` is not specified the default is `LCview`.

|Value	| Description|
|---:	|:---|
|LCView|	Default save format, does not have to be explicitly specified. Saves Q-tensor and potentials into binary files, viewable using LCView program. Files written are named as:“result<iteration number>.dat” and “result_final.dat”|
|LCviewTXT|	Saves director values and potential into human readable text files, also viewable using ||LCView. The ordering of the nodes is specified by the correponding mesh file that is also saved in the same directory. Files are named as: “result_t<iteration number>.dat” and “result_t_final.dat”.|
|RegularVTK|	Saves director order parameter and potential values interpolated onto a regular grid in a human readable text file.  Details about this format can be easily found e.g. by googling “vtk file formats”|
|RegularVecMat|	Saves director, order parameter and potential values interpolated onto a regular grid into text files that can be executed as scripts in MATLAB to load the values into variables.
|DirStackZ|	Saves director components only to a comma separated values (CSV) text file. The first line gives the number of points in x,y, and z directions as well as a time step value (in seconds).  The remaining rows list the director x,y and z components. Each row corresponds to a stack of directors along the Z-axis. Within each stack, director components are interleaved (i.e. ordered nx, ny, nz, nx, ny, nz,...). Stack locations in the x-y plane increase in the x-direction first, then y direction is increment by one.|
|VTKUnstructuredAsciiGrid|	Saves the potential values and director vector using VTK unstructured legacy grid files (ASCII text files), which also include the FE mesh at each time step. These files are comaptible for example with the ParaView data visualisation application.  For details about this format see e.g. https://kitware.github.io/vtk-examples/site/VTKFileFormats/ (“Simple Legacy Formats” section)|


### **RegularGridSize = [nX, nY, nZ]** ###
Vector specifying the size (number of points in the x, y and z directions) of a regular, evenly spaced grid used for saving results. This must be specified, for example, when `SaveFormat` is set to or includes `RegularVTK` or `DirstackZ` (or any other format that outputs to a regular grid).

The grid points are spaced evenly between the minimum and maximum coordinates in each dimension. Example:
RegularGridSize = [10,1,10] results in a 10 x 1 x 10 (2D) slice centred along mean y-coordinate of the structure.

There may be some problems with non-cuboidal structures, when regular grid points are located outside the modelling window (i.e. not in any tetrahedron). This should be a TODO item...   

### **SaveIter = 0** ###
Numerical integer value controlling the frequency, in number of iterations, of saving intermediate result files. Values larger or equal to 0 are currently supported. When `SaveIter` is 0, only initial LC configuration and final result are saved. Otherwise results are saved if the modulus of current iteration and `SaveIter` equals zero. Default value is 0.

### **SaveTime = 0.0** ###
Optional numerical value controlling the frequency in seconds of saving intermediate results.
For example, `SaveTime = 1e-6` , results in a result file being saved every 1e-6 seconds. The value is expected to be larger or equal than 0. A value of 0 means that this option is ignored. The default value is 0.

### **OutputEnergy = 0** ###
Optional numerical value determining whether LC free energy is calculated and saved in a separate matlab compatible text file (“energy.m”). 0 = no energy calculation, 1 = energy is calculated. The default value is = 0

### **StretchVector = [1.0, 1.0, 1.0]** ###
Optional vector of length 3 specifying a three dimensional scaling of the finite element mesh. All mesh node coordinates are multiplied by their corresponding `StretchVector` components. This option can be used to scale and stretch mesh files (although stretching a mesh in one direction more than other may reduce the mesh quality). The default value is `[1.0, 1.0, 1.0]`. 

### **LoadQ = result_final.dat** ###
Optional string variable specifying the filename of a previous result. The Q-tensor distribution is loaded from this file and used as a starting condition for a new simulation. The default is that no existing result is loaded.

---
## Liquid Crystal Material Parameters ##
Liquid Crystal Material parameters are defined using the key/value pairs specified  below:

Splay, twist and bend elastic coefficients:

```
	K11 = 1.0000e-11
	K22 = 1.0000e-11
	K33 = 1.0000e-11
```

Chiral pitch length, measured in metres:

```
	p0  = 0.0e-6
```

This is a optional parameter. If not assigned (or assigned 0), no chirality is assumed.

Thermotropic coefficients:

```
	A   = -1.2000e+05
	B   = -2.1333e+06
	C   = 1.7333e+06
```

Relative dielectric permittivities parallel and perpendicular to the director:

```
	eps_par = 18.5000
	eps_per = 7.0000
```
	
Flexoelectric splay and bend coefficients (which sign convention?):
	
```
	e1     = 0.0000e+00
	e3     = 0.0000e+00
```

These are optional coefficients. If not assigned (or assigned 0), no flexoelectricity is assumed.

Rotational viscosity:

```
	gamma1  = 0.0777
```

---
## Initial LC Orientation ##
Initial liquid crystal orientation is set using a number of `BOX` structures ranging from `BOX1` to `BOX99`. Each box represents a separate cuboid 3D volume whose position and size need to be specified. If box volumes are overlapping, the box with a higher number overrides the conditions of a lower numbered box within the overlapping region.

Use of Boxes is optional. When no Boxes are defined, or in regions not covered by any Boxes, a uniform LC orientation along the x-axis (no tilt, no twist) is used instead. 

**`BOXn.Type`**

Type is a string variable that defines the variation of LC within the box. Currently Type can be `Normal`, `Random`, `Hedgehog`.

* `Normal` set the director in a uniform/linearly varying orientation (variation only a function of the z-direction) 
* `Random` sets a randomized director orientation within the box.
* `Hedgehog` creates a +1 hedgehog defect at the centre of the box.

A `Type` must be defined for the n'th BOX in order for the box to be valid. If no `Type` is defined, the rest of the parameters for the n'th BOX are ignored.

**`BOXn.X, BOXn.Y, BOXn.Z`**

The size and location of each box is defined by the vectors `X`, `Y` and `Z` (of length 2, each) that range from min to max range. It is OK if the x, y, z ranges are larger than the modelled structure; the non-overlapping regions are ignored. 

**`BOXn.Tilt, BOXn.Twist`** (only used when Type is Normal)

`Tilt` and `Twist` are vector variables of length 2. `Tilt[1]` and `Twist[1]` contain the tilt and twist angles at `Z[1]`, whereas `Tilt[2]` and `Twist[2]` contain the total variations of tilt and twist within the box in the z-direction . The default values of both `Tilt` and `Twist` are `[0, 0]`, i.e. uniform LC director along the x-axis.

For example:
```
BOX1.Type = Normal
BOX1.X = [0.0, 1.0]        # a box from x = 0 to 1 
BOX1.Y = [0.0, 2.0]        #            y = 0 to 2
BOX1.Z = [0.0, 5.0]        #            z = 0 to 5
BOX1.Tilt  = [2.0, 0.0]    # tilt angle is uniform 2 degrees
BOX1.Twist = [45.0, 90.0]  # twist varies from 45  
                           # to 134 degrees along the z-axis
```
*Example of specifying an initial LC orientation region. A box of Normal type, ranging from x=0, y=0, z=0 to x=1.0, y=2.0, z=5.0 microns is defined. The tilt angle is set to 2 degrees throughout the whole of the box. The tilt angle is 45 degrees at the bottom of the box (z=0) and changes linearly by 90 degrees throughout the box, so that at the top (z=5), it is135 degrees.*

---
## Anchoring ##
Liquid crystal alignment or anchoring is set using a number of `FIXLC` structures ranging from `FIXLC1` to `FIXLC99`, as defined in the mesh.

**`FIXLCn.Anchoring`**

Anchoring is a string that describes the type of anchoring. Possible values are Homeotropic, Strong, Weak, Degenerate, Freeze or Polymerise.

|`FIXLCn.Anchoring`| Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|---:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|Strong		| Fixes the LC orientation to the easy direction.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|Weak		| Weak anchoring (strength determined by the `FIXLCn.Strength` parameter) in the easy direction.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|Homeotropic| Fixes the LC along the local surface normal.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|Degenerate | The LC is free to rotate in the local surface plane. If `FIXLCn.Strength` is less than 0, this becomes weak Homeotropic anchoring.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|Freeze		| Fixes the LC to its initial orientation. Can be used e.g. to create non-uniform alignment.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|Polymerise	| Special undocumented secret feature that does not require an alignment surface in the mesh! Fixes the value of the Q-tensor at any nodes where the order parameter value is equal or below the value specified in `FIXLCn.Strength.`                                                                                                                                                                                                                                                                                                                                                                                         |
|ManualNodes| Use the `FIXLCn.Params` vector to manually specify mesh node numbers where to strongly fix the LC orientation to a given easy direction. You may have to look up the node numbering in the mesh file. **Note**: mesh node numbering in `qlc3d` starts at 0 whereas in mesh file it may start at 1, so subtract one, and in meshes with dielectric volume regions, node numbering is reordered in qlc3d so that all LC nodes come before dielectric nodes. This means that it is safer to look up node numbering from a previously written result mesh file than from one written by GiD or Gmsh when dielectrics are present.|

**`FIXLCn.Strength`**

`Strength` is a numeric variable that defines the Rapini-Papoular Anchoring strength when Weak or Degenerate anchoring is used. If Strength is set to a negative value and Type is Degenerate, a weak homeotropic anchoring with strength |`Strength`| is achieved.

**`FIXLCn.Easy`**

`Easy` is a numeric vector of length 3. It sets the easy tilt, twist and rotation (around the director) angles of the anchoring. The rotation angle is ignored for all cases except when `Anchoring` is `Weak` and `FIXLCn.K1` is not equal to `FIXLCn.K2`.

**`FIXLCn.K1, FIXLCn.K2`**

`K1` and `K2` are optional numeric values that can be used with `Anchoring` = `Weak` to specify anisotropic anchoring strengths. The actual polar and azimuthal anchoring strengths are `Strength * K1 `and `Strength * K2` respectively. Default values are 1 for both, resulting in isotropic anchoring, where there is no difference in the in-plane vs. out-of-plane anchoring strengths.

**`FIXLCn.Params = [ list of numbers ]`**

A list of additional options can be provided to different anchoring types. The meaning of the values will depend on the corresponding anchoring type.

When `FIXLCn.Anchoring = ManualNodes`, the `Params` are mesh node indexes. In other cases this is ignored and not needed. 

**`FIXLCn.overrideVolume = true | false`**

A boolean flag that controls setting of the initial LC orientation at the FixLC surface. When set to `true`, the LC initial orientation
at the surface is set to the surface's easy direction. When set to `false`, the initial LC orientation at the surface is set
according to the settings for the volume in that region, e.g. using `Box`es or by loading an initial configuration from an
existing result file using `LoadQ` parameter.

By default, this is set to `true`, so that it should only be necessary to use this if you want to change it to `false`. 

### Example
```
FIXLC1.Anchoring = Homeotropic
FIXLC1.Strength = 1e-4
FIXLC1.Easy = [80.0000, 45.0000, 0.0000]
FIXLC1.K1 = 1.0000
FIXLC1.K2 = 1.0000
```
*Example of specifying anchoring conditions. The anchoring on FixLC1 surface is set to strong homeotropic (perpendicular to the surface plane). In this case, the parameters Strength, Easy , K1 and K2 are ignored, since homeotropic anchoring type overrides these settings.*

---
## Electric Potentials ##
The potentials on electrodes are defined using a number of data structures ranging from `E1` to `E99`. Each structure corresponds to an electrode as defined in the mesh.

**`En.Time`**

`Time` is a numerical vector. Each value in `Time[i]` corresponds to a switching of potential `Pot[i]` on the electrodes. **Note**, the length of `Time` should be equal to the length of `Pot`.

**`En.Pot`**

`Pot` is a numerical vector. Each value in `Pot[i]` corresponds to a potential switching at time `Time[i]`. **Note**, the length of `Pot` should be equal to the length of `Time`.

### Example ###
```
E1.Time = [0.0, 1e-3]
E1.Pot  = [3.0, 0.0]
```
*Example of defining electrode potentials. Electrode 1 is switched to 3 V at time 0ms and to 0V at 1ms.*

---
## Uniform Electric Fields ##
A uniform electric field may be defined by setting the x, y and z-components of the EField vector. The units of the vector are V/μm. For example, the following defines a constant field of 1.3 V / μm throughout the structure.

```
EField = [0 , 0 , 1.3]
```

Uniform electric fields can be defined in structures that do not contain any electrodes. Uniform electric fields are not functions of time, i.e. swithing on/off during a simulation is not currently supported.

What should happen when both electrodes and a uniform E-Field have been defined? Clearly we can't have both at the same time!

---
## Dielectric Properties ##

eps_dielectric is a numerical vector that holds the relative permittivity values of dielectric regions defined in the structure.

For example, the permittinvity in a structure containing two dielectric regions, `Dielectric1` and `Dielectric2`, is defined as 
```
eps_dielectric = [1.0,5.0]
```
---
## Mesh Adaptation ##

Mesh adaptation is optional, and is controlled by specifying `REFINEMENT` objects (see section **REFINEMENT Objects** below), and the instances in time and/or iterations (see section **Periodically Repeating Mesh Adaptation**) when mesh adaptation occurs.

### Periodically Repeating Mesh Adaptation ###
The frequency of repeated mesh adaptation is specified by assigning a value to either variable `RefRefIter` or `RepRefTime`:
```
	RepRefIter = iteration number
```
or
```
	RepRefTime = time period [milliseconds]. # not working!!
```

If `RepRefIter` or `RepRefTime` are not specified, or if they are both set to a zero value, the mesh is not adapted periodically. However, it is still possible to specify explicit iteration/time instances to do this in the individual REFINEMENT objects.

Note that if a REFINEMENT object specifies iteration numbers or instances in time, these will take precedence over RepRefIter/RepRefTime for the corresponding REFINMENT object.

### REFINEMENT Objects ###

`REFINEMENT` objects are used to specify different properties (method, region, etc.) of mesh adaptation. Currently up to 99 `REFINEMENT` objects are supported. The following format is used:
```
REFINEMENTi.<key> = <value>
```
Where i is a number in the range 1-99, specifying the `REFINEMENT` object number, <key> is a property name string and <value> is the corresponding value (string, number or array).

The following key-value pairs can be defined:

|Key|Value Type|Description|
|---:|:---:|:---|
|Type|String|Sets the type of the corresponding REFINEMENT object. This will define how the other key/value pairs are interpreted for this object. See examples below for possible values.|
|Values|Array of numbers|Interpretation depends on Type. See Examples below.|
|Iterations|Array of numbers|Explicit Iteration numbers when adaptation described by this object is performed.|
|Times|Array of numbers|Explicit times when adaptation described by this object is performed.|
|X|Array of numbers|Array of x-coordinates|
|Y|Array of numbers|Array of y-coordinates|
|Z|Array of numbers|Array of z-coordinates|


### `REFINEMENTi.Type = Sphere` ###

Selects elements within a spherical region centred on X, Y and Z. The Values array contains the radius/radii of the region(s). 

Example:
```
	REFINEMENT1.Type = Sphere
	REFINEMENT1.Values = [0.5, 0.25, 0.1]
	REFINEMENT1.X = [0.5]
	REFINEMENT1.Y = [0.5]
	REFINEMENT1.Z = [0.5]
```
Performs three refinement iterations, selecting elements at 0.5, 0.25 and 0.1 micron distances from coordinate 0.5, 0.5, 0.5 (i.e. hree concentric spherical regions are selected).

### `REFINEMENTi.Type=Box` ###
Selects elements within cuboidal regions specified in the X, Y and Z arrays.

Example:
```
	REFINEMENT1.Type = Box
	REFINEMENT1.X = [0, 0.1, 0, 0.05]
	REFINEMENT1.Y = [0, 0.1, 0, 0.05]
	REFINEMENT1.Z = [0, 0.1, 0, 0.05]
```
Performs two refinement iterations, selecting elements within the two cubes ranging from x, y, z = 0 to 0.1 microns for the first cube and x,y,z = 0 to 0.05 microns for the second one (the second cube is inside the first one in this case).

### `REFINEMENTi.Type = Change` ###

Selects elements in which any of the five Q-tensor components change is equal or larger than specified in the Values array.

Example:
```
	REFINEMENT1.Type = Change
	REFINEMENT1.Values = [0.2, 0.2]
```
Performs two refinement iterations, selecting any elements within the whole modelled structure where the Q-tensor change exceeds the value 0.2. This refinement is applied globally to any elements in the whole structure.

---
## Solver Settings ##
Various settings can be used to controll some numerical and implementational aspects of a simulation. In most cases there should be no need to set or change these values, the `qlc3d` application should be able to choose appropriate settings automatically. However, advanced users or developers may want to override some of these settings. 

The equations for the Q-tensor and the electric potential are solved using either the Preconditioned Conjugate Gradient (**PCG**) or the Gerneralised Minimal Residual (**GMRES**) method. In most cases, the **PCG** algorithm is faster and requires less memory, but it requires the solution matrix to be symmetric, which is not the case for certain combinations of liquid crystal material parameters. To manually choose between these algorithms, set `Q_Solver` and `V_Solver` to either 1 or 0. If the algorithm is not explicitly defined, `qlc3d` automatically chooses the appropriate one depending on whether the solution matrix is symmetric.

Both the **PCG** and **GMRES** algorithms require a preconditioner matrix that can be diagonal (Jacobi), incomplete cholesky or incomplete LU decomposition. Se the Q and V preconditioner variables to 1, 2, or 3. 

A maximum number of iterations and required accuracy must also be defined for both methods.

Usually the **PCG** method with diagonal preconditioning is fastest, but in some cases does not converge to a solution. The **GMRES** method with e.g. incomplete LU preconditioning is more robust, but a bit slower and uses more memory. For example, structures with Neumann boundary conditions require using the **GMRES** method for the potential solution since the matrix problem is non-symmetric. (detection of cases like this should really be automatic…)

The dynamics of the Q-tensor involves multiple sub-iterations per time step. If the size of the time step is too large, these sub-iterations may not converge. The `Q_Newton_Panic_Iter` determines the maximum number of sub iterations before the time step is scaled by the value of `Q_Newton_Panic_Coeff`.

```
	nThreads 	= 4
	Q_Solver	= 1	# 0 = PCG, 1 = GMRES
	V_Solver	= 1	# 0 = PCG, 1 = GMRES
#--------------------------------------------------------
#	Q-Tensor solver settings
#--------------------------------------------------------
	Q_Newton_Panic_Iter  = 10
	Q_Newton_Panic_Coeff = 0.1		
	
	# Preconditioned Conjugate Gradient Q-tensor solver settings
	Q_PCG_Preconditioner 	= 0	# 0 = Diagonal
						# 1 = Incomplete Cholesky
						# 2 = Incoplete LU
	Q_PCG_Maxiter	= 2000
	Q_PCG_Toler		= 1e-3
		
	# GMRES Q-tensor solver settings	
	Q_GMRES_Preconditioner	= 2	# 0 = Diagonal
							# 1 = Incomplete Cholesky
							# 2 = Incoplete LU
	Q_GMRES_Maxiter	= 200
	Q_GMRES_Restart	= 100
	Q_GMRES_Toler	= 1e-3
#-------------------------------------------------------
#	Potential solver settings
#-------------------------------------------------------
	# Preconditioned Conjugate Gradient potential solver settings
	V_PCG_Preconditioner = 0		# 0 = Diagonal
							# 1 = Incomplete Cholesky
							# 2 = Incomplete LU
	V_PCG_Maxiter	= 2000
	V_PCG_Toler		= 1e-6
	
	# GMRES potential solver settings	
	V_GMRES_Preconditioner	= 2		# 0 = Diagonal
							# 1 = Incomplete Cholesky
							# 2 = Incoplete LU
	V_GMRES_Maxiter	= 2000
	V_GMRES_Restart	= 50;
	V_GMRES_Toler	= 1e-6;
```

---
## Settings File Example ##
See the examples included in this project for more examples of settings files. Below, a "standalone" settings file is given that can be used as a starting point in a project.


```
meshName = mesh3df.txt 

# Simulate 10 milliseconds
EndCriterion = Time
EndValue = 10e-3

dt = 1e-9 #dt>0 for time stepping, dt=0 for Newton method
dtLimits = [1e-9, 1e-4]   # min, max allowed dt
MaxError = 1e-5 # time step accuracy
TargetdQ = 1e-2 # predictor-corrector param.
DtFunction = [0.25, 0.8, 1.2, 10]; 

Savedir = results_directory # results are saved here
SaveIter = 5 # write result to file every 5 iterations

K11 = 6.2e-12 # LC elastic coefficients
K22 = 3.9e-12
K33 = 8.2e-12
p0 =  2.0e-6  # chiral pitch length is 2 microns

A = -0.240000e6 # -2 degrees
B = -2.133300e6
C = 1.733300e6
eps_par = 18.5; # LC dielectric coefficients
eps_per = 7.0;
gamma1 = 0.0777; # LC rotational viscosity

# set initial LC orientation to 10 degrees tilt everywhere
bOX1.Type = Normal
BOX1.Params = [1,0]
BOX1.X = [0.0, 50]
BOX1.Y = [0.0, 50]
BOX1.Z = [0, 1]
BOX1.Tilt = [10,0]
BOX1.Twist = [0,0]

# FIXLC1 surface at 10 degrees tilt, strong anchoring
FIXLC1.anchoring =strong 
FIXLC1.Strength = 1e-4
FIXLC1.K1 = 1.0
FIXLC1.K2 = 1.0
FIXLC1.Easy = [10, 0, 0]

# FIXLC2 surface at 1 degree tilt, 90 degrees twist
FIXLC2.anchoring =strong 
FIXLC2.Strength = 1e-4
FIXLC2.K1 = 1.0
FIXLC2.K2 = 1.0
FIXLC2.Easy = [1.0, 90, 0]

# specify potentials for two electrodes
E1.Time = [0]
E1.Pot = [0]

E2.Time = [0, 5e-3] # 5V at start, but switch off to 0 at 5 milliseconds
E2.Pot = [5, 0]

# mesh contains 2 different dielectric regions. 
# Set their permittivities to 4 and 5
eps_dielectric = [4, 5];

# Do mesh refinement every 30 iterations
RepRefIter = 30
#split elements where Q changes by more than 0.1 per element.
REFINEMENT1.Type = Change
REFINEMENT1.Values = [0.1, 0.1] # Q change threshold 
REFINEMENT1.X = [0, 1, 0, 0.2]
REFINEMENT1.Y = [0, 1, 0, 1]
REFINEMENT1.Z = [0, 10, 0, 0.5]
```

---
## Minimal Mesh File Example ##
For additional information about meshing and geometry, see [this](./mesh.md).

Below a minimal GiD mesh file is shown. It's not useful for any realistic simulations, but is provided here as an example of the expected Gid mesh file format. 

If your own GiD mesh file looks like gibberish when opened in a text editor, it may have been saved in binary format and qlc3d will not be able to read it. To save in text format in GiD, a mesh must be *exported*, not *saved*.

* The below mesh consists of 9 nodes in range 0 to 1 microns (section “Coordinates”). 
* It has 12 tetrahedral elements (first section “Elements”), all of LC material (material number 4). 
* It has periodic surfaces consisting of 6 x prism elements (second section “Elements”).
* It has 12 triangular elements (third section “Elements”) of material number 3 (periodic). 
* The mesh does not have any alignment surfaces nor electrodes, i.e. it is a fully periodic box of LC material.

```
MESH dimension 3 ElemType Tetrahedra Nnode 4
Coordinates
    1               1               0               1
    2     0.499998865     0.500003276     0.499995589
    3               1               1               1
    4               0               0               1
    5               1               0               0
    6               1               1               0
    7               0               0               0
    8               0               1               1
    9               0               1               0
end coordinates

Elements
       1       9       8       4       2 4
       2       4       3       1       2 4
       3       4       8       3       2 4
       4       5       1       6       2 4
       5       7       4       1       2 4
       6       9       4       7       2 4
       7       6       9       7       2 4
       8       5       6       7       2 4
       9       7       1       5       2 4
      10       1       3       6       2 4
      11       3       8       9       2 4
      12       6       3       9       2 4
END ELEMENTS
MESH dimension 3 ElemType Prism Nnode 6
Coordinates
end coordinates

Elements
      13       8       4       9       3       1       6
      14       4       7       9       1       5       6
      15       3       1       4       6       5       7
      16       8       3       4       9       6       7
      17       4       1       7       8       3       9
      18       1       5       7       3       6       9
end elements
MESH dimension 3 ElemType Triangle Nnode 3
Coordinates
end coordinates

Elements
      19       8       4       9 3
      20       4       7       9 3
      21       4       1       7 3
      22       1       5       7 3
      23       3       6       1 3
      24       1       6       5 3
      25       8       9       3 3
      26       3       9       6 3
      27       3       1       4 3
      28       8       3       4 3
      29       6       7       5 3
      30       9       7       6 3
end elements
```