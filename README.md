# README #

### TO DO LIST PHYSICS ###

1. Add elasticity
1. Add plasticity
1. Add self gravity
1. Compute Shear Heating and add it to temperature
1. Compute Viscous flow law according to the viscous strain rate, as opposed to the total strain rate
1. Add compressibility

### TO DO LIST OPTIMIZATION ###

1. optimize marker to cells to make less boolean tests
1. add subgrid diffusion on markers
1. Physics_InterpFromParticlesToCell: values of T might need to be overwritten for periodic nodes
1. Create Numbering->boundary and Numbering->boundaryNeighbour to optimize the looping the boundaries for EC arrays
1. For periodic boundaries instead of numbering lines like this: 0 1 2 3 0 1, change to numebring: 3 0 1 2 3 1, or when computing interpolation to Boundaries use the solution from the inner node
1. Remove particles from the air

### TO DO LIST VISUALIZATION ###
hCalling Python to do data treatment would be ideal. See :ttp://www.linuxjournal.com/article/8497

### TO DO LIST SETUP ###
1. Add sandbox type BC including constant heat flux boundaries
1. Read geomIO file and perform point in polygon test

### DONE ###
- Stokes solver
- Non linear rheology
- Penalty method
- Pure shear BC
- Periodic BC
- Linked list of markers
- Velocity, strain rate and viscosity visualization
- Interpolation of Viscosity from markers to cell centers
- Interpolation of Viscosity froms cell centers to markers
- Interpolation routine from cell centers to nodes
- Add gravity
- Use texture instead of triangles
- Visualize markers using a geometry shader or point sprites
- Add passive grid
- Add temperature