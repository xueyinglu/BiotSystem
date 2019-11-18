#include "BiotSystem.h"

#include "PressureSolution.h"
#include "DisplacementSolution.h"
using namespace std;
void BiotSystem::calc_error()
{
    Vector<float> difference_per_cell_pressure(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler_pressure,
                                      solution_pressure,
                                      PressureSolution(t),
                                      difference_per_cell_pressure,
                                      QGauss<dim>(fe_pressure.degree + 2),
                                      VectorTools::L2_norm);
    double L2_norm_pressure = difference_per_cell_pressure.l2_norm();
        //VectorTools::compute_global_error(triangulation,
        //                                  difference_per_cell_pressure,
        //                                  VectorTools::L2_norm);
    /* Calculate the L2 norm of displacement solution */

    Vector<float> difference_per_cell_displacement(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler_displacement,
                                      solution_displacement,
                                      DisplacementSolution(t),
                                      difference_per_cell_displacement,
                                      QGauss<dim>(fe_displacement.degree + 2),
                                      VectorTools::L2_norm);
    double L2_norm_displacement = difference_per_cell_displacement.l2_norm();

    /*
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();
    const FEValuesExtractors::Vector displacements (0);
    const unsigned int n_q_points = quadrature_displacement.size();
    vector<Vector<double>> sol_u_values(n_q_points,Vector<double> (dim));
    vector<Tensor<1,dim>> ana_u_values(n_q_points);
    DisplacementSolution disp(t);

    double L2_norm_displacement = 0;
    for(; cell!= endc; ++cell){
        fe_value_displacement.reinit(cell);
        fe_value_displacement.get_function_values(solution_displacement, sol_u_values);
        disp.value_list(fe_value_displacement.get_quadrature_points(), ana_u_values);

        for (unsigned int q =0; q<n_q_points; q++){
            Tensor<1, dim> sol;
            sol[0] =sol_u_values[q][0];
            sol[1] = sol_u_values[q][1];
            Tensor<1,dim> cell_difference= sol - ana_u_values[q];
            L2_norm_displacement += cell_difference.norm_square() *fe_value_displacement.JxW(q);
        }
    }
    */
    L2_norm_displacement = sqrt(L2_norm_displacement);

    convergence_table.add_value("time", t);
    convergence_table.add_value("1/h", 1./h);
    convergence_table.add_value("L2_u", L2_norm_displacement);
    convergence_table.add_value("L2_p", L2_norm_pressure);
}