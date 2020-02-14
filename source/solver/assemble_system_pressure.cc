#include "BiotSystem.h"
#include "AuxTools.h"
#include "TerzaghiPressure.h"
#include "PermFunction.h"
using namespace std;
void BiotSystem::assemble_system_pressure()
{
    system_matrix_pressure.reinit(sparse_pattern_pressure);
    system_rhs_pressure.reinit(dof_handler_pressure.n_dofs());
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
    PermFunction perm_function;
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
        if (test_case == TestCase::heterogeneous){
            perm_function.value_list(fe_value.get_quadrature_points(), permeability_values);
        }
        fe_value.get_function_values(prev_timestep_sol_pressure, prev_timestep_sol_pressure_values);
        fe_value.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        fe_value_displacement.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_sol_grad_u_values);
        fe_value_displacement.get_function_gradients(prev_fs_sol_displacement, prev_fs_sol_grad_u_values);
        /*
        cout << "perm =" << permeability_values[0] << endl;
        cout << "mu_f =" << mu_f << endl;
        cout << "biot_alpha = " << biot_alpha << endl;
        cout << "biot_inv_M =" << biot_inv_M <<endl;
        cout << "lame_lambda = " << lame_lambda <<endl;
        cout << "lame_mu = " << lame_mu <<endl;
        cout << "K_b = " << K_b << endl;
        */
        /* assemble cell level matrix and rhs */
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            // calculate the mean stress values at the quadrature point
            const Tensor<2, dim> prev_time_grad_u = Tensors ::get_grad_u<dim>(q, prev_timestep_sol_grad_u_values);
            const double prev_time_div_u = Tensors ::get_divergence_u<dim>(prev_time_grad_u);
            const Tensor<2, dim> prev_fs_grad_u = Tensors ::get_grad_u<dim>(q, prev_fs_sol_grad_u_values);
            const double prev_fs_div_u = Tensors ::get_divergence_u<dim>(prev_fs_grad_u);
            prev_timestep_mean_stress = K_b * prev_time_div_u - biot_alpha * prev_timestep_sol_pressure_values[q];
            prev_fs_mean_stress = K_b * prev_fs_div_u - biot_alpha * prev_fs_sol_pressure_values[q];
            for (unsigned int i = 0; i < dofs_per_cell; i++)
            {
                for (unsigned int j = 0; j < dofs_per_cell; j++)
                {
                    // elliptic part
                    cell_matrix(i, j) +=
                        (del_t / mu_f * permeability_values[q] * // 1/mu_f * k
                         fe_value.shape_grad(i, q) *             // grad phi_i(x_q)
                         fe_value.shape_grad(j, q) *             // grad phi_j(x_q)
                         fe_value.JxW(q));                       // dx
                    // parabolic part
                    cell_matrix(i, j) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) *                              // (1/M + alpha^2/K_b)/del_t
                         fe_value.shape_value(i, q) * fe_value.shape_value(j, q) * fe_value.JxW(q)); // phi(x_q)*phi(x_q) dx
                }
                // source term
                //cell_rhs(i) +=
                //    (fe_value.shape_value(i, q) * // phi_i(x_q)
                //     1 *                          // f(x_q)
                //     fe_value.JxW(q));            // dx

                // prev time step
                cell_rhs(i) +=
                    ((biot_inv_M + biot_alpha * biot_alpha / K_b) * // (1/M + alpha^2/K_b)/del_t
                     prev_timestep_sol_pressure_values[q] *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));

                // change in mean stress
                cell_rhs(i) -=
                    (biot_alpha / K_b *
                     (prev_fs_mean_stress - prev_timestep_mean_stress) *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        constraints_pressure.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix_pressure);
        constraints_pressure.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs_pressure);
        /*
        for (int i = 0; i < dofs_per_cell; i++)
        {
            for (int j = 0; j < dofs_per_cell; j++)
            {
                system_matrix_pressure.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
            }
            system_rhs_pressure(local_dof_indices[i]) += cell_rhs(i);
        }
        */
    }
    if (test_case == TestCase::benchmark)
    {
        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 0,
                                                 ZeroFunction<dim>(),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
    }
    else if (test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous)
    { 
        cout << "Terzaghi pressure BC" << endl;
        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 2, // boundary id
                                                 ConstantFunction<2>(pressure_dirichlet_bc),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
        
        /*
        cout << "------- DEBUG Terzaghi pressure --------- " << endl;
        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 0,
                                                 TerzaghiPressure(t),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 1,
                                                 TerzaghiPressure(t),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 2,
                                                 TerzaghiPressure(t),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
        */
    }
}