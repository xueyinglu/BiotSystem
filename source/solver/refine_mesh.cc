#include "BiotSystem.h"

void BiotSystem::refine_mesh()
{

    /*
    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      cell_eta_u,
                                                      0.6,
                                                      0.4);
                                                      */
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, cell_eta_p,0.1, 0.4);                                                 
    const unsigned int max_grid_level = 9;
    const unsigned int min_grid_level = 3;
    if (triangulation.n_levels() > max_grid_level)
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
    Vector<double> prev_sol_p = solution_pressure;
    Vector<double> prev_sol_u = solution_displacement;
    solution_trans_p.prepare_for_coarsening_and_refinement(prev_sol_p);
    solution_trans_u.prepare_for_coarsening_and_refinement(prev_sol_u);
    triangulation.execute_coarsening_and_refinement();
    setup_system();
    solution_trans_p.interpolate(prev_sol_p, solution_pressure);
    solution_trans_u.interpolate(prev_sol_u, solution_displacement);
    constraints_pressure.distribute(solution_pressure);
    constraints_displacement.distribute(solution_displacement);

}