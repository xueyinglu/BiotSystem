#include "BiotSystem.h"
#include "PressureSolution.h"
using namespace std;
void BiotSystem::run_fixed_stress()
{
    make_grid();
    setup_system();
    // Initialize u_0
    cout<<"Initializing u_0"<<endl;
    assemble_system_displacement();
    solve_displacement();
    prev_timestep_sol_displacement = solution_displacement;
    // interpolate initial pressure
    VectorTools::interpolate(dof_handler_pressure,
                            PressureSolution(0),
                            prev_timestep_sol_pressure);
    
    for (timestep = 1; timestep < (T / del_t); timestep++)
    // for (timestep = 1; timestep < 2 ; timestep++)
    {
        cout << "timestep = " << timestep << endl;
        t += del_t;
        fixed_stress_iteration();
        output_displacement(timestep, -1);
        output_pressure(timestep, -1);
        calc_error();
        calc_a_posteriori_indicators();
        prev_timestep_sol_displacement = solution_displacement;
        prev_timestep_sol_pressure = solution_pressure;
    }
    output_error();
}