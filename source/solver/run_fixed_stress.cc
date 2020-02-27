#include "BiotSystem.h"
#include "PressureSolution.h"
#include "TerzaghiPressure.h"
using namespace std;
void BiotSystem::run_fixed_stress()
{
    set_material_properties();
    make_grid();
    setup_system();
    // interpolate initial pressure
    if (test_case == TestCase::benchmark)
    {
        VectorTools::interpolate(dof_handler_pressure,
                                 PressureSolution(0),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        // Initialize u_0
        cout << "Initializing u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
        prev_timestep_sol_displacement = solution_displacement;
    }
    else if (test_case == TestCase::terzaghi)
    { // p_0 = 0; u_0 = 0;
        cout << "Benchmark Terzaghi : p_0 =0; u_0 = 0" << endl;
        VectorTools::interpolate(dof_handler_pressure,
                                 ZeroFunction<dim>(),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        VectorTools::interpolate(dof_handler_displacement,
                                 ZeroFunction<dim>(dim),
                                 solution_displacement);
        prev_timestep_sol_displacement = solution_displacement;
    }
    else if (test_case == TestCase::heterogeneous)
    {
        cout << "heterogeneous test" << endl;

        VectorTools::interpolate(dof_handler_pressure,
                                 ConstantFunction<dim>(initial_pressure_value),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        cout << "Initializing u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
        prev_timestep_sol_displacement = solution_displacement;
    }
    initial_pressure = solution_pressure;
    initial_displacement = solution_displacement;

    for (timestep = 1; timestep < ((T +1e-5)/ del_t); timestep++)
    {
        cout << "timestep = " << timestep << endl;
        t += del_t;
        fixed_stress_iteration();
        plot_error();
        // output_displacement(timestep, -1);
        // output_pressure(timestep, -1);
        if (test_case == TestCase::benchmark || test_case == TestCase::terzaghi)
        {
            calc_error();
        }

        if (criteria != 3)
        {
            calc_a_posteriori_indicators_p();
            calc_a_posteriori_indicators_u();
        }
        if (test_case == TestCase::benchmark)
        {
            calc_efficiency_indices();
        }
        calc_strain_stress();
        if (adaptivity == true)
        {
            refine_mesh();
        }

        prev_timestep_sol_displacement = solution_displacement;
        prev_timestep_sol_pressure = solution_pressure;
    }
    output_error();
}
