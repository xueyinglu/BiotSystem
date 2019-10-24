#include "BiotSystem.h"
#include "RightHandSide.h"
#include "InitialPressure.h"
void BiotSystem::assemble_system_displacement()
{
    QGauss<dim> quadrature_formula(2);

    FEValues<dim> fe_values(fe_displacement, quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);
    // need this to access pressure solutions
    FEValues<dim> fe_values_pressure(fe_pressure, quadrature_formula,
                                     update_values |
                                         update_quadrature_points |
                                         update_JxW_values |
                                         update_gradients);

    const unsigned int dofs_per_cell = fe_displacement.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    ConstantFunction<dim> lambda(1.), mu(1.);
    RightHandSide right_hand_side;
    InitialPressure initial_pressure;

    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);
    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(dim));
    std::vector<double> pore_pressure_values(n_q_points);

    // Now we can begin with the loop over all cells:
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_pressure = dof_handler_pressure.begin_active();

    for (; cell != endc; ++cell, ++cell_pressure)
    {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values_pressure.reinit(cell_pressure);

        // Next we get the values of the coefficients at the quadrature
        // points. Likewise for the right hand side:
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);
        if (time_step == 0)
        {
            initial_pressure.value_list(fe_values_pressure.get_quadrature_points(), pore_pressure_values);
        }
        else
        {
            fe_values_pressure.get_function_values(solution_pressure, pore_pressure_values);
        }
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            const unsigned int
                component_i = fe_displacement.system_to_component_index(i).first;

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
                const unsigned int
                    component_j = fe_displacement.system_to_component_index(j).first;

                for (unsigned int q_point = 0; q_point < n_q_points;
                     ++q_point)
                {
                    cell_matrix(i, j) +=
                        ((fe_values.shape_grad(i, q_point)[component_i] *
                          fe_values.shape_grad(j, q_point)[component_j] *
                          lambda_values[q_point]) +
                         (fe_values.shape_grad(i, q_point)[component_j] *
                          fe_values.shape_grad(j, q_point)[component_i] *
                          mu_values[q_point]) +
                         ((component_i == component_j) ? (fe_values.shape_grad(i, q_point) *
                                                          fe_values.shape_grad(j, q_point) *
                                                          mu_values[q_point])
                                                       : 0)) *
                        fe_values.JxW(q_point);
                }
            }
        }

        // Assembling the right hand side is also just as discussed in the
        // introduction:
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            const unsigned int
                component_i = fe_displacement.system_to_component_index(i).first;

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
                cell_rhs(i) += fe_values.shape_value(i, q_point) *
                               rhs_values[q_point](component_i) *
                               fe_values.JxW(q_point);
                // coupling pressure
                cell_rhs(i) += fe_values.shape_value(i, q_point) *
                               biot_alpha * pore_pressure_values[q_point] *
                               fe_values.JxW(q_point);
            }
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix_displacement.add(local_dof_indices[i],
                                               local_dof_indices[j],
                                               cell_matrix(i, j));

            system_rhs_displacement(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    hanging_node_constraints.condense(system_matrix_displacement);
    hanging_node_constraints.condense(system_rhs_displacement);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             0,
                                             ZeroFunction<dim>(dim),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
}