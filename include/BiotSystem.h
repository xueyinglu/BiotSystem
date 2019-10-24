#ifndef BIOT_SYSTEM_H_
#define BIOT_SYSTEM_H_
#include "DealiiHeader.h"

using namespace dealii;

class BiotSystem
{
public:
    BiotSystem();
    // virtual BiotSystem();
    void run_fixed_stress();

private:
    double del_t = 0.1;
    double T = 1;
    int time_step = 0;

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

    // coupling

    double biot_alpha = 0.8;
    double fs_threshold = 1e-8;
    Vector<double> prev_timestep_sol_pressure;
    Vector<double> prev_timestep_sol_displacement;
    Vector<double> prev_fs_sol_pressure;
    Vector<double> prev_fs_sol_displacement;

    void make_grid();
    void setup_system();

    void assemble_system_pressure();
    void assemble_system_displacement();

    void solve_pressure();
    void solve_displacement();

    void fixed_stress_iteration();

    void output_results() const;

};

#endif