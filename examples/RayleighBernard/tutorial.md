# This tutorial contains the Rayleigh-Bernard convection in a rectangular domain

## Step 1: Generating the mesh using the gmsh script

```
gmsh RB.geo -3 -o RB.msh
```

## Step 2: Decomposing the mesh using the Morton (z-order) space filling curve

a. Make a copy of the input file and name it decompose.input (name does not matter and can be different)
b. Single rank runs dont require step 2a
c. Replace the preprocessed field `Block A` with `Block B`

Block A
```
&dopamine_domain
  mesh_file = 'mesh'
  mesh_type = 'preprocessed'
/
```

Block B
```
&dopamine_domain
  mesh_file = 'RB.msh'
  mesh_type = 'gmsh'
/
```

d. Decompose the mesh using the dopamine-decompose utility

```
dopamine-decompose decompose.input --nranks 8 --output-dir mesh
```

## Step 3: Running the solver

```
mpirun -np 8 dopamine RB.input > dope_run.log 2>&1 &
disown
```

## Step 4: Visualisation

Based on the write frequency of 3D snapshots, the solver should generate a `RB.xmf` file which can be imported into paraview to visualise the results.

Animation created using the helper python script (`RBAnalysis.py`)

![test](../../figures/rb.gif)