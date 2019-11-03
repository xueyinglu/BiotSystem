#include "BiotSystem.h"
#include "PressureSolution.h"
void BiotSystem::plot_error() const{
    Vector<double> interpolated_exact_sol(dof_handler_pressure.n_dofs());
    Vector<double> error(dof_handler_pressure.n_dofs());
    VectorTools::interpolate (dof_handler_pressure,
                              PressureSolution(t),
                              interpolated_exact_sol);
    error = interpolated_exact_sol;
    error -= solution_pressure;
    for (int i = 0; i < error.size(); i++){
        error[i] = std::abs(error[i]);
    }
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_pressure);
    data_out.add_data_vector(solution_pressure, "pressure");
    data_out.add_data_vector(interpolated_exact_sol, "exact_sol");
    data_out.add_data_vector(error, "error_p");
    data_out.build_patches();
    std::ofstream output("output/error-" + std::to_string(timestep) +".vtk");
    data_out.write_vtk(output);
}

