# Disc quadrilateral mesh
Implementation of the algorithm propused in "A PDE based approach to multi-domain
partitioning and quadrilateral meshing" (see [here](https://www.ljll.math.upmc.fr/~frey/publications/imr21.pdf)). Although the implementation is general
it has only been tested on a semicircle. The inputs are the radius of the circle, the center and the number of subdivisions of the centre square.

# Result
The following is the main result of the implementation. The resulting mesh can be used for any purpose. 
<img src="mesh.jpg" width="500">
# Pedagogical value
I think the main contribution of the project is its pedagogical value. Through the code you can inspect all the main parts of the complex algorithm.
These parts are shown below. 
## 1. Creation of the boundary vector field
<img src="images/boundary_vf.png" width="500">

## 2. Mesh grid and delaunay triangulation
<img src="images/delaunay_triangulation.png" width="500">

## 3. FEM solution
<img src="images/FEM_solution.png" width="500">

## 4. Find the singular points and the cross field
<img src="images/singular_pts.png" width="500">

## 5. Find the streamlines by advection
<img src="images/advection.png" width="500">

## 6. The final mesh is obtained using the transfinite algorithm in the streamlines
<img src="images/transfinite.png" width="500">
