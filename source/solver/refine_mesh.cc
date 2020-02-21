#include "BiotSystem.h"

void BiotSystem::refine_mesh()
{

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      cell_eta_u,
                                                      0.6,
                                                      0.4);
    const unsigned int max_grid_level = 9;
    const unsigned int min_grid_level = 3;
    if (triangulation.n_levels > max_grid_level)
    {
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(max_grid_level))
            cell->clear_refine_flag();
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(min_grid_level))
            cell->clear_coarsen_flag();
    }
    SolutionTransfer<dim> solution_trans_p(dof_handler_pressure);
    SolutionTransfer<dim> solution_trans_u(dof_handler_displacement);
    triangulation.prepare_coarsening_and_refinement();
    solution_trans_p.prepare_coarsening_and_refinement(prev_timestep_sol_pressure);
    solution_trans_p.prepare_coarsening_and_refinement(prev_timestep_sol_displacement);
}