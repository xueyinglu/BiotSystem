#include "BiotSystem.h"
#include "PressureSolution.h"
#include "DisplacementSolution.h"
#include <iostream>
#include <fstream>
using namespace std;
void BiotSystem::check_disp_solver_convergence()
{
    num_global_refinement = 1;

    GridGenerator::hyper_cube(triangulation, 0, 1);
    cout << "================================displacement solver convergence test================================" << endl;
    for (int cycle = 1; cycle < 8; cycle++)
    {
        cout << "cycle = " << cycle << endl;
        num_global_refinement *= 2;
        h = 1. / num_global_refinement;
        // make_grid();
        triangulation.refine_global(1);
        setup_system();
        VectorTools::interpolate(dof_handler_pressure,
                                 PressureSolution(0),
                                 solution_pressure);
        assemble_system_displacement();
        solve_displacement();
        QGauss<dim> quadrature_displacement(fe_displacement.degree + 1);
        FEValues<dim> fe_value_displacement(fe_displacement,
                                            quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                       endc = dof_handler_displacement.end();
        const unsigned int n_q_points = quadrature_displacement.size();
        vector<Vector<double>> sol_u_values(n_q_points, Vector<double>(dim));
        vector<Tensor<1, dim>> ana_u_values(n_q_points);
        DisplacementSolution disp(0);

        double L2_norm_displacement = 0;
        for (; cell != endc; ++cell)
        {
            fe_value_displacement.reinit(cell);
            fe_value_displacement.get_function_values(solution_displacement, sol_u_values);
            disp.value_list(fe_value_displacement.get_quadrature_points(), ana_u_values);

            for (unsigned int q = 0; q < n_q_points; q++)
            {
                Tensor<1, dim> sol;
                sol[0] = sol_u_values[q][0];
                sol[1] = sol_u_values[q][1];
                Tensor<1, dim> cell_difference = sol - ana_u_values[q];
                L2_norm_displacement += cell_difference.norm_square()* fe_value_displacement.JxW(q);
            }
        }
        L2_norm_displacement = sqrt(L2_norm_displacement);
        cout << "L2_norm_displacement = " << L2_norm_displacement << endl;

        convergence_table.add_value("cycle", cycle);
        convergence_table.add_value("1/h", 1. / h);
        convergence_table.add_value("L2_u", L2_norm_displacement);
    }
    ofstream u_convergence_file("u_convergence.tex");
    convergence_table.write_tex(u_convergence_file);
}