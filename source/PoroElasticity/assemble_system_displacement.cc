#include "BiotSystem.h"
#include "DisplacementSolution.h"
#include "AuxTools.h"
using namespace std;
void BiotSystem::assemble_system_displacement()
{
    QGauss<dim> quadrature_formula(fe_displacement.degree + 1);

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

    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);
    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(dim));
    //std::vector<Vector<double>> grad_p_values(n_q_points,
    //                                       Vector<double>(dim));
    std::vector<double> pore_pressure_values(n_q_points);
    std::vector<Tensor<1, dim>> grad_p_values(n_q_points);
    Tensor<2,dim> identity = Tensors::get_Identity<dim> ();
    const FEValuesExtractors::Vector displacements (0);
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
        /*
        if (timestep == 0) // initialize u_0
        {
            initial_pressure.value_list(fe_values_pressure.get_quadrature_points(), pore_pressure_values);
            //initial_pressure.grad_list(fe_values_pressure.get_quadrature_points(), grad_p_values);
            // for (int q = 0 ; q <n_q_points;q++){
            //     cout << "initial pressure [q] = " <<pore_pressure_values[q]<<endl;
            //}
        }
        else
        {
            fe_values_pressure.get_function_values(solution_pressure, pore_pressure_values);
        }
        */

        fe_values_pressure.get_function_values(solution_pressure, pore_pressure_values);
        fe_values_pressure.get_function_gradients(solution_pressure, grad_p_values);
        //initial_pressure.value_list(fe_values_pressure.get_quadrature_points(), pore_pressure_values);
       /* 
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
        
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            const unsigned int
                component_i = fe_displacement.system_to_component_index(i).first;

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            { // body force
                // cell_rhs(i) += fe_values.shape_value(i, q_point) *
                //                rhs_values[q_point](component_i) *
                //               fe_values.JxW(q_point);
                // coupling pressure
                cell_rhs(i) += fe_values.shape_grad(i, q_point)[component_i] *
                               biot_alpha * pore_pressure_values[q_point] *
                               fe_values.JxW(q_point);
                //cell_rhs(i) -= fe_values.shape_value(i, q_point) *
                //               biot_alpha * grad_p_values[q_point][component_i] *
                 //              fe_values.JxW(q_point);
            }
        }
        */
        
        // Assemble the cell matrix as in elasticity_cg
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            std::vector<Tensor<1, dim>> phi_i_u(dofs_per_cell);
            std::vector<Tensor<2, dim>> phi_i_grads_u(dofs_per_cell);
            std::vector<Tensor<2, dim>> E_phi(dofs_per_cell);
            std::vector<Tensor<2, dim>> sigma_phi(dofs_per_cell);

            // Compute and store desired quantities
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                phi_i_u[k] =fe_values[displacements].value(k,q);
                phi_i_grads_u[k] = fe_values[displacements].gradient(k, q);
                E_phi[k] = 0.5 * (phi_i_grads_u[k] + transpose(phi_i_grads_u[k]));
                sigma_phi[k] = 2.0 * mu_values[q] * E_phi[k] + lambda_values[q] * trace(E_phi[k]) * identity;
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                { // Displacements
                    cell_matrix(j, i) += fe_values.JxW(q) * scalar_product(sigma_phi[i], E_phi[j]);
                }

                // assemble cell level rhs as in elasticity_cg
                 cell_rhs(i) += biot_alpha * pore_pressure_values[q] * trace(phi_i_grads_u[i]) *fe_values.JxW(q);
                // cell_rhs(i) -= biot_alpha * (grad_p_values[q]*phi_i_u[i])*fe_values.JxW(q);
            }
            
        } // end q_point
        
        
        cell->get_dof_indices(local_dof_indices);
        /*
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix_displacement.add(local_dof_indices[i],
                                               local_dof_indices[j],
                                               cell_matrix(i, j));

            system_rhs_displacement(local_dof_indices[i]) += cell_rhs(i);
        }
        */
        hanging_node_constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix_displacement);
        hanging_node_constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs_displacement);
    }

    // hanging_node_constraints.condense(system_matrix_displacement);
    // hanging_node_constraints.condense(system_rhs_displacement);
    /*
    vector<bool> component_mask;
    component_mask.push_back(false);
    component_mask.push_back(true);
    std::map<types::global_dof_index, double> boundary_values; */
    /*VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             0,
                                             ZeroFunction<dim>(dim),
                                             boundary_values); */
    /*
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             1,
                                             ZeroFunction<dim>(dim),
                                             boundary_values,
                                             ComponentMask(component_mask));

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
    component_mask[0] = true;
    component_mask[1] = false;
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             0,
                                             ZeroFunction<dim>(dim),
                                             boundary_values,
                                             ComponentMask(component_mask));
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
    */
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             0,
                                             DisplacementSolution(t),
                                             //ZeroFunction<dim>(dim),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
    
}