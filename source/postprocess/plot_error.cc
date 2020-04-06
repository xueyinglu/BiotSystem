#include "BiotSystem.h"
#include "PressureSolution.h"
#include "TerzaghiPressure.h"
#include "DisplacementSolution.h"
using namespace std;
void BiotSystem::plot_error() const
{
    if (test_case == benchmark || test_case == terzaghi)
    {
        Vector<double> interpolated_exact_sol(dof_handler_pressure.n_dofs());
        Vector<double> error(dof_handler_pressure.n_dofs());
        if (test_case == TestCase::benchmark)
        {
            VectorTools::interpolate(dof_handler_pressure,
                                     PressureSolution(t),
                                     interpolated_exact_sol);
        }
        else if (test_case == TestCase::terzaghi)
        {
            VectorTools::interpolate(dof_handler_pressure,
                                     TerzaghiPressure(t),
                                     interpolated_exact_sol);
        }

        error = interpolated_exact_sol;
        error -= solution_pressure;
        for (int i = 0; i < error.size(); i++)
        {
            error[i] = std::abs(error[i]);
        }
        Vector<double> interpolated_exact_sol_u(dof_handler_displacement.n_dofs());
        Vector<double> error_u(dof_handler_displacement.n_dofs());
        VectorTools::interpolate(dof_handler_displacement,
                                 DisplacementSolution(t),
                                 interpolated_exact_sol_u);
        error_u = interpolated_exact_sol_u;
        error_u -= solution_displacement;
        for (int i = 0; i < error_u.size(); i++)
        {
            error_u[i] = std::abs(error_u[i]);
        }
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler_pressure);
        data_out.add_data_vector(solution_pressure, "pressure");
        data_out.add_data_vector(interpolated_exact_sol, "exact_sol");
        data_out.add_data_vector(error, "error_p");
        data_out.build_patches();
        ofstream output("visual/error-p-" + std::to_string(timestep) + ".vtk");
        data_out.write_vtk(output);

        DataOut<dim> data_out_u;
        data_out_u.attach_dof_handler(dof_handler_displacement);
        vector<string> u_error_names;
        u_error_names.push_back("u_x_error");
        u_error_names.push_back("u_y_error");
        vector<string> u_names;
        u_names.push_back("u_x");
        u_names.push_back("u_y");
        vector<string> u_exact_names;
        u_exact_names.push_back("u_x_exact");
        u_exact_names.push_back("u_y_exact");
        data_out_u.add_data_vector(error_u, u_error_names);
        data_out_u.add_data_vector(solution_displacement, u_names);
        data_out_u.add_data_vector(interpolated_exact_sol_u, u_exact_names);
        data_out_u.build_patches();
        ofstream output_u("visual/error-u-" + std::to_string(timestep) + ".vtk");
        data_out_u.write_vtk(output_u);
    }

    else if (test_case == TestCase::heterogeneous)
    {
        DataOut<dim> data_out_prop;
        data_out_prop.attach_dof_handler(dof_handler_output);
        Vector<double> cell_perm;
        cell_perm.reinit(dof_handler_output.n_dofs());
        VectorTools::interpolate(dof_handler_output,
                                 perm_function,
                                 cell_perm);
        Vector<double> cell_lambda;
        cell_lambda.reinit(dof_handler_output.n_dofs());
        VectorTools::interpolate(dof_handler_output,
                                 lambda_function,
                                 cell_lambda);
        Vector<double> cell_mu;
        cell_mu.reinit(dof_handler_output.n_dofs());
        VectorTools::interpolate(dof_handler_output,
                                 mu_function,
                                 cell_mu);

        //data_out_prop.add_data_vector(cell_perm, "perm", DataOut<dim>::type_dof_data);
        //data_out_prop.add_data_vector(cell_lambda, "lambda", DataOut<dim>::type_dof_data);
        //data_out_prop.add_data_vector(cell_mu, "mu", DataOut<dim>::type_dof_data);
        cout << "line 93" << endl;
        data_out_prop.add_data_vector(cell_perm, "perm",DataOut<dim>::type_dof_data);
        data_out_prop.add_data_vector(cell_lambda, "lambda",DataOut<dim>::type_dof_data);
        data_out_prop.add_data_vector(cell_mu, "mu",DataOut<dim>::type_dof_data);
        data_out_prop.build_patches();
        ofstream output_prop("visual/property-" + std::to_string(timestep) + ".vtk");
        data_out_prop.write_vtk(output_prop);

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler_pressure);
        data_out.add_data_vector(solution_pressure, "p");
        data_out.build_patches();
        ofstream output("visual/" + filename_base + "-p-" + std::to_string(timestep) + ".vtk");
        data_out.write_vtk(output);

        DataOut<dim> data_out_u;
        data_out_u.attach_dof_handler(dof_handler_displacement);
        vector<string> u_names;
        u_names.push_back("u_x");
        u_names.push_back("u_y");
        data_out_u.add_data_vector(solution_displacement, u_names);
        data_out_u.build_patches();
        ofstream output_u("visual/" + filename_base + "-u-" + std::to_string(timestep) + ".vtk");
        data_out_u.write_vtk(output_u);
    }
}
