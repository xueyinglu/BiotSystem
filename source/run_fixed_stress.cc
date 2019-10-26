#include "BiotSystem.h"
using namespace std;
void BiotSystem::run_fixed_stress(){
    make_grid();
    setup_system();
    // Initialize u_0
    assemble_system_displacement();
    solve_displacement();
    prev_timestep_sol_displacement = solution_displacement;

    for (timestep = 0 ; timestep < (T/del_t ) ; timestep++){
        cout << "timestep = " <<timestep << endl;
        fixed_stress_iteration();
        prev_timestep_sol_displacement = solution_displacement;
        prev_timestep_sol_pressure = solution_pressure;
    }
    output_results();
}