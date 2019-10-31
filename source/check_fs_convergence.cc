#include "BiotSystem.h"
using namespace std;
// check the convergence of fixed stress, return the residual
double BiotSystem::check_fs_convergence(int fs_count)
{
    QGauss<dim> quadrature_pressure(fe_pressure.degree + 2);
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_pressure(fe_pressure,
                                    quadrature_pressure, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    const unsigned int n_q_points = quadrature_pressure.size();

    vector<vector<Tensor<1, dim>>> prev_fs_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<double> prev_fs_sol_pressure_values(n_q_points);
    vector<double> pressure_values(n_q_points);

    double mean_stress;
    double prev_fs_mean_stress;
    double residual = 0.0;

    for (; cell != endc; ++cell, ++cell_displacement)
    {
        fe_value_pressure.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        fe_value_pressure.get_function_values(solution_pressure, pressure_values);
        if (timestep == 1 && fs_count == 1)
        {
            initial_pressure.value_list(fe_value_pressure.get_quadrature_points(), prev_fs_sol_pressure_values);
        }
        else
        {
            fe_value_pressure.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        }
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        fe_value_displacement.get_function_gradients(prev_fs_sol_displacement, prev_fs_sol_grad_u_values);

        for (unsigned int q = 0; q < n_q_points; q++)
        {

            mean_stress = K_b * (grad_u_values[q][0][0] + grad_u_values[q][1][1]) - biot_alpha * pressure_values[q];
            prev_fs_mean_stress = K_b * (prev_fs_sol_grad_u_values[q][0][0] + prev_fs_sol_grad_u_values[q][1][1]) - biot_alpha * prev_fs_sol_pressure_values[q];
            residual += (mean_stress - prev_fs_mean_stress) * (mean_stress - prev_fs_mean_stress) * fe_value_pressure.JxW(q);
        }
    }
    cout << "fixed stress iteration convergence criteria = " << sqrt(residual) << endl;
    return sqrt(residual);
}