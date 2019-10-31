#include "BiotSystem.h"

void BiotSystem::make_grid(){
    GridGenerator::hyper_cube(triangulation, 0, 1);
    //triangulation.begin_active() -> face(0)-> set_boundary_id(1); 
    // mark the 0-face of the coarsest element with boundary_id =1. Neumann BC will be applied to this face.
    triangulation.refine_global(num_global_refinement);
    std::cout <<"Number of active cells:" << triangulation.n_active_cells() << std::endl;
}