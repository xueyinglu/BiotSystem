#include "BiotSystem.h"
using namespace std;
void BiotSystem::calc_a_posteriori_indicators()
{
    /* calculate eta_alg */
    double eta_alg_m;
    double dum1 = 0;
    double dum2 = 0;
    for (vector<double>::iterator it = eta_fs.begin(); it != eta_fs.end(); it++)
    {
        dum1 += *it;
        dum2 += ((*it) * (*it));
    }
    eta_alg_m = del_t * dum1 * dum1 + h * h * dum2;
    eta_alg.push_back(eta_alg_m);

    /* calculate eta_time */
    double eta_t_p_n = 0; //eta_t_p_n^2
    QGauss<dim> quadrature_pressure(fe_pressure.degree + 2);
    FEValues<dim> fe_value_pressure(fe_pressure, quadrature_pressure,
                                    update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    const unsigned int n_q_points = quadrature_pressure.size();
    vector<Tensor<1, dim>> grad_p_values(n_q_points);
    vector<Tensor<1, dim>> prev_timestep_grad_p_values(n_q_points);
    vector<double> permeability_values(n_q_points);
    for (; cell != endc; ++cell)
    {
        fe_value_pressure.reinit(cell);
        fe_value_pressure.get_function_gradients(solution_pressure, grad_p_values);
        if (timestep == 1)
        {
            initial_pressure.grad_list(fe_value_pressure.get_quadrature_points(), prev_timestep_grad_p_values);
        }
        else
        {
            fe_value_pressure.get_function_gradients(prev_timestep_sol_pressure, prev_timestep_grad_p_values);
        }
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<1, dim> cell_difference = grad_p_values[q] - prev_timestep_grad_p_values[q];
            eta_t_p_n += permeability_values[q] *
                         cell_difference.norm_square() *
                         fe_value_pressure.JxW(q);
        }
    }
    eta_t_p_n = eta_t_p_n * del_t / 3;
    cout << "eta_t_p_n squared = " << eta_t_p_n << endl;
    double eta_time_m;
    if (timestep == 1)
    {
        eta_time_m = eta_t_p_n;
        eta_time.push_back(eta_time_m);
    }
    else
    {
        eta_time_m = eta_time.back() + eta_t_p_n;
        eta_time.push_back(eta_time_m);
    }

    /* calculate eta_flow */

    double eta_flow_m = 0;
    // residual at cells
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement, quadrature_displacement,
                                        update_values | update_quadrature_points | update_gradients | update_JxW_values);
    vector<vector<Tensor<1, dim>>> prev_timestep_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));

    vector<double> prev_timestep_p_values(n_q_points);
    vector<double> p_values(n_q_points);
    vector<double> laplacian_p_values(n_q_points);
    cell = dof_handler_pressure.begin_active();
    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    for (; cell != endc; ++cell, ++cell_displacement)
    {
        fe_value_pressure.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        fe_value_pressure.get_function_values(solution_pressure, p_values);
        fe_value_pressure.get_function_values(prev_timestep_sol_pressure, prev_timestep_p_values);
        fe_value_pressure.get_function_laplacians(solution_pressure, laplacian_p_values);
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        fe_value_displacement.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_grad_u_values);
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        for (unsigned int q = 0; q < n_q_points; q++)
        {   
            // TODO: extension to non constant permeability
            double cell_residual = 1 / mu_f * permeability_values[q] * laplacian_p_values[q] 
                                - biot_inv_M / del_t * (p_values[q] - prev_timestep_p_values[q]) 
                                - biot_alpha / del_t * (grad_u_values[q][0][0] + grad_u_values[q][1][1] 
                                - prev_timestep_grad_u_values[q][0][0] - prev_timestep_grad_u_values[q][1][1]);
            // cout << "cell residual = " << cell_residual << endl;
            eta_flow_m += cell_residual * cell_residual * fe_value_pressure.JxW(q);
        }
    }

    eta_flow_m = sqrt(eta_flow_m) * h * h * del_t;
    /*
    // jump of flux at cell boundaries.
    double flux_jump = 0;
    QGauss<dim> face_quadrature(fe_pressure.degree + 2);
    const unsigned int n_face_q_points = face_quadrature.size();
    FEFaceValues<dim> fe_face_values(fe_pressure, face_quadrature,
                                    update_values | update_normal_vectors
                                    | update_gradients | update_quadrature_points
                                    | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_values(fe_pressure, face_quadrature,
                                    update_values | update_normal_vectors
                                    | update_gradients | update_quadrature_points
                                    | update_JxW_values);
    
    
    vector<Tensor<1, dim>> face_grad_p_values(n_face_q_points);
    cell = dof_handler_pressure.begin_active();
    for (; cell != endc; ++cell){
        fe_value_pressure.reinit(cell);
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no){
            
            
            const unsigned int neighbor_face = cell -> neighbor_of_neighbor(face_no);
            fe_face_values.reinit(cell, face_no);
            fe_face_neighbor_values.re
            fe_face_values.get_function_gradients(solution_pressure, face_grad_p_values);
            for (unsigned int q = 0; q < n_face_q_points; q++){
                flux_jump += 
            }
        }
    }
    */
    a_posterior_indicators_table.add_value("time", t);
    a_posterior_indicators_table.add_value("eta_fs", eta_fs.back());
    a_posterior_indicators_table.add_value("eta_alg", eta_alg_m);
    a_posterior_indicators_table.add_value("eta_time", eta_time_m);
    a_posterior_indicators_table.add_value("eta_flow", eta_flow_m);
}