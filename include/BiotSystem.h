#ifndef BIOT_SYSTEM_H_
#define BIOT_SYSTEM_H_
#include "DealiiHeader.h"
#include "InitialPressure.h"
#include "RightHandSide.h"
using namespace dealii;
using namespace std;
class BiotSystem
{
public:
    BiotSystem();
    // virtual BiotSystem();
    void run_fixed_stress();

private:
    double del_t = 2.5e-4;
    double T = 0.05;
    double t = 0;
    int timestep = 0;
    int num_global_refinement = 2;
    double h = 1./num_global_refinement;

    Triangulation<dim> triangulation;
    // pressure solution
    FE_Q<dim> fe_pressure;
    DoFHandler<dim> dof_handler_pressure;
    SparsityPattern sparse_pattern_pressure;
    SparseMatrix<double> system_matrix_pressure;

    Vector<double> solution_pressure;
    Vector<double> system_rhs_pressure;

    // displacement solution
    FESystem<dim> fe_displacement;
    ConstraintMatrix hanging_node_constraints;
    DoFHandler<dim> dof_handler_displacement;

    SparsityPattern sparsity_pattern_displacement;
    SparseMatrix<double> system_matrix_displacement;

    Vector<double> solution_displacement;
    Vector<double> system_rhs_displacement;
    ConvergenceTable convergence_table;

    // Data
    double mu_f = 1; // fluid viscosity
    RightHandSide right_hand_side; // mechanics equation body force
    InitialPressure initial_pressure;
    ConstantFunction<dim> permeability;
    ConstantFunction<dim> lambda, mu;
    // coupling

    double biot_alpha = 0.75;
    double K_b = 0.5;
    double biot_inv_M = 3./28;
    double tol_fixed_stress = 1e-5;
    Vector<double> prev_timestep_sol_pressure;
    Vector<double> prev_timestep_sol_displacement;
    Vector<double> prev_fs_sol_pressure;
    Vector<double> prev_fs_sol_displacement;

    /* global a posteriori error estimators (recorded for each time step) */
    vector<double> eta_fs;
    vector<double> eta_alg;
    vector<double> eta_time;
    vector<double> eta_flow;

    ConvergenceTable a_posterior_indicators_table;


    void make_grid();
    void setup_system();

    void assemble_system_pressure(int fs_count);
    void assemble_system_displacement();

    void solve_pressure();
    void solve_displacement();

    void fixed_stress_iteration();

    double check_fs_convergence(int fs_count); // check the convergence of fixed-stress iteration

    void output_displacement(int timestep, int fs_count) const;
    void output_pressure(int timestep, int fs_count) const;
    void output_error();
    void calc_error(); // compute the errors
    void process_solution(int fs_count); // compute the errors
    void plot_error(int fs_count) const;

    void calc_a_posteriori_indicators();

};

#endif