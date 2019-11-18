#include "BiotSystem.h"
using namespace std;
void BiotSystem::assemble_system_pressure()
{
    QGauss<dim> quadrature(fe_pressure.degree + 1);
    FEValues<dim> fe_value(fe_pressure,
                           quadrature, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    const unsigned int dofs_per_cell = fe_pressure.dofs_per_cell;
    const unsigned int n_q_points = quadrature.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    vector<double> permeability_values(n_q_points);
    vector<double> prev_timestep_sol_pressure_values(n_q_points);
    vector<double> prev_fs_sol_pressure_values(n_q_points);
    vector<vector<Tensor<1, dim>>> prev_timestep_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> prev_fs_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    double prev_timestep_mean_stress;
    double prev_fs_mean_stress;
    // vector<vector<double>> grad_u(dim, dim);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    for (; cell != endc; ++cell, ++cell_displacement)
    {
        fe_value.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        cell_matrix = 0;
        cell_rhs = 0;
        /* get the function values at current element */
        permeability.value_list(fe_value.get_quadrature_points(), permeability_values);
        /*
        if (timestep == 1)
        {
            initial_pressure.value_list(fe_value.get_quadrature_points(), prev_timestep_sol_pressure_values);
            if (fs_count == 1)
            {
                initial_pressure.value_list(fe_value.get_quadrature_points(), prev_fs_sol_pressure_values);
            }
            else
            {

                fe_value.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
            }
        }
        else
        {
            fe_value.get_function_values(prev_timestep_sol_pressure, prev_timestep_sol_pressure_values);
            fe_value.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        }
        */
        fe_value.get_function_values(prev_timestep_sol_pressure, prev_timestep_sol_pressure_values);
        fe_value.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        fe_value_displacement.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_sol_grad_u_values);
        fe_value_displacement.get_function_gradients(prev_fs_sol_displacement, prev_fs_sol_grad_u_values);
        /* assemble cell level matrix and rhs */
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            // calculate the mean stress value at the quadrature point
            prev_timestep_mean_stress = K_b * (prev_timestep_sol_grad_u_values[q][0][0] + prev_timestep_sol_grad_u_values[q][1][1])
                                         - biot_alpha * prev_timestep_sol_pressure_values[q];
            prev_fs_mean_stress = K_b * (prev_fs_sol_grad_u_values[q][0][0] + prev_fs_sol_grad_u_values[q][1][1])
                                         - biot_alpha * prev_fs_sol_pressure_values[q];
            for (unsigned int i = 0; i < dofs_per_cell; i++)
            {
                for (unsigned int j = 0; j < dofs_per_cell; j++)
                {
                    // elliptic part
                    cell_matrix(i, j) +=
                        (1. / mu_f * permeability_values[q] * // 1/mu_f * k
                         fe_value.shape_grad(i, q) *          // grad phi_i(x_q)
                         fe_value.shape_grad(j, q) *          // grad phi_j(x_q)
                         fe_value.JxW(q));                    // dx
                    // parabolic part
                    cell_matrix(i, j) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t *                      // (1/M + alpha^2/K_b)/del_t
                         fe_value.shape_value(i, q) * fe_value.shape_value(j, q) * fe_value.JxW(q)); // phi(x_q)*phi(x_q) dx
                }
                // source term
                //cell_rhs(i) +=
                //    (fe_value.shape_value(i, q) * // phi_i(x_q)
                //     1 *                          // f(x_q)
                //     fe_value.JxW(q));            // dx

                // prev time step
                cell_rhs(i) +=
                    ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t *                      // (1/M + alpha^2/K_b)/del_t
                     prev_timestep_sol_pressure_values[q] *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));

                // change in mean stress
                cell_rhs(i) -=
                    (biot_alpha / K_b / del_t *
                     (prev_fs_mean_stress - prev_timestep_mean_stress) *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; i++)
        {
            for (unsigned int j = 0; j < dofs_per_cell; j++)
            {
                system_matrix_pressure.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
            }
            system_rhs_pressure(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    std::map<types::global_dof_index, double> boundary_values;
    //VectorTools:: interpolate_boundary_values(dof_handler_pressure, 0, Functions::ZeroFunction<dim>(), boundary_values);
    VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                             0,
                                             ZeroFunction<2>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
}