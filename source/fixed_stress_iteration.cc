#include "BiotSystem.h"
using namespace std;
void BiotSystem::fixed_stress_iteration()
{
    bool iteration = true;
    int fs_count = 0;
    double l2_norm;
    prev_fs_sol_displacement = prev_timestep_sol_displacement;
    prev_fs_sol_pressure = prev_timestep_sol_pressure;
    while (iteration)
    {
        fs_count++;
        cout << "fixed stress no. " << fs_count << endl;
        assemble_system_pressure();
        solve_pressure();
        assemble_system_displacement();
        solve_displacement();
        l2_norm = check_fs_convergence(fs_count);
        iteration = (l2_norm > tol_fixed_stress);

        prev_fs_sol_pressure = solution_pressure;
        prev_fs_sol_displacement = solution_displacement;
        // process_solution(fs_count);
        // plot_error(fs_count);
    }

    eta_fs.push_back(1 / sqrt(del_t) * l2_norm);
}