#include "BiotSystem.h"
using namespace std;
void BiotSystem::fixed_stress_iteration()
{
    bool iteration = true;
    int fs_count = 0;
    double l2_norm;
    while (iteration)
    {
        prev_fs_sol_pressure = solution_pressure;
        prev_fs_sol_displacement = solution_displacement;
        fs_count++;
        cout << "fixed stress no. " << fs_count << endl;

        assemble_system_pressure();

        solve_pressure();
        assemble_system_displacement();
        solve_displacement();
        if (criteria !=3)
        {
            l2_norm = check_fs_convergence();
            iteration = (l2_norm > tol_fixed_stress);
        }
        else
        {
            l2_norm = check_fs_convergence();
            eta_fs.push_back(1 / sqrt(del_t) * l2_norm); 
            calc_a_posteriori_indicators_p();
            calc_a_posteriori_indicators_u();
            iteration = (eta_alg.back() > 0.1 * (eta_flow.back() + eta_time.back() 
            + eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back() + eta_u.back()));
            if (iteration){
                eta_fs.pop_back();
                eta_alg.pop_back();
                eta_flow.pop_back();
                eta_time.pop_back();
                eta_face_partial_sigma_n.pop_back();
                eta_N_partial_sigma_n.pop_back();
                eta_face_partial_sigma.pop_back();
                eta_face_sigma.pop_back();
                eta_face_sigma_n.pop_back();
                eta_N_sigma_n.pop_back();
                eta_u.pop_back();
                eta_u_n.pop_back();

            }
        }

        // process_solution(fs_count);
        // plot_error(fs_count);
    }
    num_fs.push_back(fs_count);
    //TODO: If convergence criteria = 2 this is not correct
    eta_fs.push_back(1 / sqrt(del_t) * l2_norm);
}