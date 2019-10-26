#include "BiotSystem.h"
using namespace std;
void BiotSystem::fixed_stress_iteration(){
    bool iteration = true;
    int fs_count = 0;
    prev_fs_sol_displacement = prev_timestep_sol_displacement;
    prev_fs_sol_pressure = prev_timestep_sol_pressure;
    while(iteration){
        fs_count ++;
        cout<<"fixed stress no. " << fs_count <<endl;
        assemble_system_pressure(fs_count);
        solve_pressure();
        assemble_system_displacement();
        solve_displacement();
        iteration = convergence_fixed_stress(fs_count);
        prev_fs_sol_pressure = solution_pressure;
        prev_fs_sol_displacement = solution_displacement;
    }
}