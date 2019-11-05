#include "BiotSystem.h"
#include "AuxTools.h"
using namespace std;
#include <numeric>
void BiotSystem::calc_a_posteriori_indicators_u()
{

    /* calculate eta_face_partial_sigma_n */
    double eta_face_partial_sigma_t = 0;
    QGauss<dim - 1> face_quadrature(fe_pressure.degree + 1);
    // const unsigned int n_face_q_points = face_quadrature.size();
    FEFaceValues<dim> fe_face_p(fe_pressure, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_p(fe_pressure, face_quadrature,
                                         update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_u(fe_displacement, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_u(fe_displacement, face_quadrature,
                                         update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe_displacement.dofs_per_cell;
    Tensor<2, dim> identity = Tensors::get_Identity<dim>();
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    typename DoFHandler<dim>::active_cell_iterator
        cell_u = dof_handler_displacement.begin_active();

    for (; cell != endc; ++cell, ++cell_u)
    {
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            const auto face_p = cell->neighbor(face_no);
            const auto face_u = cell_u->neighbor(face_no);
            if (!face_p->at_boundary())
            {
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                const auto neighbor_p = cell->neighbor(face_no);
                const auto neighbor_u = cell_u->neighbor(face_no);
                vector<vector<Tensor<1, dim>>> face_grad_u_values(fe_face_u.n_quadrature_points);
                vector<vector<Tensor<1, dim>>> neighbor_grad_u_values(fe_face_u.n_quadrature_points);
                vector<vector<Tensor<1, dim>>> prev_timestep_face_grad_u_values(fe_face_u.n_quadrature_points);
                vector<vector<Tensor<1, dim>>> prev_timestep_neighbor_grad_u_values(fe_face_u.n_quadrature_points);
                vector<double> face_p_values(fe_face_p.n_quadrature_points);
                vector<double> neighbor_p_values(fe_face_p.n_quadrature_points);
                vector<double> prev_timestep_face_p_values(fe_face_p.n_quadrature_points);
                vector<double> prev_timestep_neighbor_p_values(fe_face_p.n_quadrature_points);
                vector<double> lambda_values(fe_face_p.n_quadrature_points);
                vector<double> mu_values(fe_face_p.n_quadrature_points);

                const unsigned int neighbor_face_p = cell->neighbor_of_neighbor(face_no);
                const unsigned int neighbor_face_u = cell_u->neighbor_of_neighbor(face_no);
                fe_face_p.reinit(cell, face_no);
                fe_face_neighbor_p.reinit(neighbor_p, neighbor_face_p);
                fe_face_u.reinit(cell_u, face_no);
                fe_face_neighbor_u.reinit(neighbor_u, neighbor_face_u);

                fe_face_p.get_function_values(solution_pressure, face_p_values);
                fe_face_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_face_p_values);
                fe_face_neighbor_p.get_function_values(solution_pressure, neighbor_p_values);
                fe_face_neighbor_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_neighbor_p_values);

                fe_face_u.get_function_gradients(solution_displacement, face_grad_u_values);
                fe_face_u.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_face_grad_u_values);
                fe_face_neighbor_u.get_function_gradients(solution_pressure, neighbor_grad_u_values);
                fe_face_neighbor_u.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_neighbor_grad_u_values);

                lambda.value_list(fe_face_p.get_quadrature_points(), lambda_values);
                mu.value_list(fe_face_p.get_quadrature_points(), mu_values);
                vector<Point<dim>> v_normal1 = fe_face_p.get_normal_vectors();
                vector<Point<dim>> v_normal2 = fe_face_neighbor_p.get_normal_vectors();

                for (unsigned int q = 0; q < fe_face_p.n_quadrature_points; q++)
                {
                    Tensor<2, dim> face_grad_u = Tensors::get_grad_u<dim>(q, face_grad_u_values) - Tensors::get_grad_u<dim>(q, prev_timestep_face_grad_u_values);
                    Tensor<2, dim> face_E = 0.5 * (face_grad_u + transpose(face_grad_u));
                    Tensor<2, dim> face_sigma = 2 * mu_values[q] * face_E + lambda_values[q] * trace(face_E) * identity;

                    Tensor<2, dim> neighbor_grad_u = Tensors::get_grad_u<dim>(q, neighbor_grad_u_values) - Tensors::get_grad_u<dim>(q, prev_timestep_neighbor_grad_u_values);
                    Tensor<2, dim> neighbor_E = 0.5 * (neighbor_grad_u + transpose(neighbor_grad_u));
                    Tensor<2, dim> neighbor_sigma = 2 * mu_values[q] * neighbor_E + lambda_values[q] * trace(neighbor_E) * identity;

                    Tensor<1, dim> normal1, normal2;
                    normal1[0] = v_normal1[q](0);
                    normal1[1] = v_normal1[q](1);
                    normal2[0] = v_normal2[q](0);
                    normal2[1] = v_normal2[q](1);

                    Tensor<1, dim> dum = (face_sigma + biot_alpha * (face_p_values[q] - prev_timestep_face_p_values[q]) * identity) * normal1 
                                        - (neighbor_sigma + biot_alpha * (neighbor_p_values[q] - prev_timestep_neighbor_p_values[q]) * identity) * normal2;
                    eta_face_partial_sigma_t += dum.norm_square() * fe_face_p.JxW(q);
                }
            }
        }
    }

    eta_face_partial_sigma_t = sqrt(h * eta_face_partial_sigma_t);
    if (timestep == 1)
    {
        eta_face_sigma.push_back(eta_face_partial_sigma_t * eta_face_partial_sigma_t);
    }
    else
    {
        eta_face_sigma.push_back(eta_face_partial_sigma_n.back() * eta_face_partial_sigma_n.back() + eta_face_partial_sigma_t * eta_face_partial_sigma_t);
    }
    eta_face_partial_sigma_n.push_back(eta_face_partial_sigma_t);
    double dum2 = 0;
    for (auto &n : eta_face_partial_sigma_n)
    {
        dum2 += n;
    }
    eta_face_partial_sigma.push_back(dum2 * dum2);
}