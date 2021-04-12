# MATLAB binary generation
gen_mesh is top level call -- makes use of JdR .nc file for bounds. Set grid and tile arrangement here.
calls:
1) readDataLarge_pahol.m: interpolates output to horizontal grid (need to set paths)
2) rdmds_init.m: interpolates to vertical grid (give grid size here)
3) gendata.m: initialises shelf topography, shelf mass, bed and initial eta (optional)
