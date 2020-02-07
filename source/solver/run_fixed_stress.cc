#include "BiotSystem.h"
#include "PressureSolution.h"
using namespace std;
void BiotSystem::run_fixed_stress()
{
    make_grid();
    setup_system();
    // interpolate initial pressure

    VectorTools::interpolate(dof_handler_pressure,
                             PressureSolution(0),
                             solution_pressure);
    prev_timestep_sol_pressure = solution_pressure;
    // Initialize u_0
    cout << "Initializing u_0" << endl;
    assemble_system_displacement();
    solve_displacement();
    prev_timestep_sol_displacement = solution_displacement;
    for (timestep = 1; timestep < (T / del_t); timestep++)
    {
        cout << "timestep = " << timestep << endl;
        t += del_t;
        fixed_stress_iteration();
        // output_displacement(timestep, -1);
        // output_pressure(timestep, -1);
        plot_error();
        calc_error();
        
        if (criteria != 3){
           calc_a_posteriori_indicators_p();
           calc_a_posteriori_indicators_u();
        }
        
        calc_efficiency_indices();
        prev_timestep_sol_displacement = solution_displacement;
        prev_timestep_sol_pressure = solution_pressure;
    }
    output_error();
}