#include "BiotSystem.h"

void BiotSystem::setup_system()
{
    dof_handler_pressure.distribute_dofs(fe_pressure);
    std::cout << "Number of degrees of freedom for pressure: " << dof_handler_pressure.n_dofs() << std::endl;

    DynamicSparsityPattern dsp_pressure(dof_handler_pressure.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp_pressure);
    sparse_pattern_pressure.copy_from(dsp_pressure);

    system_matrix_pressure.reinit(sparse_pattern_pressure);
    solution_pressure.reinit(dof_handler_pressure.n_dofs());
    system_rhs_pressure.reinit(dof_handler_pressure.n_dofs());

    dof_handler_displacement.distribute_dofs(fe_displacement);
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler_displacement,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    DynamicSparsityPattern dsp_displacement(dof_handler_displacement.n_dofs(), dof_handler_displacement.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_displacement, dsp_displacement,
                                    hanging_node_constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern_displacement.copy_from(dsp_displacement);

    system_matrix_displacement.reinit(sparsity_pattern_displacement);

    solution_displacement.reinit(dof_handler_displacement.n_dofs());
    system_rhs_displacement.reinit(dof_handler_displacement.n_dofs());
}