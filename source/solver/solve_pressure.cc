#include "BiotSystem.h"
void BiotSystem::solve_pressure(){
    SolverControl solver_control(1000, 1e-12);
    SolverCG<> solver(solver_control);
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix_pressure, 1.2);
    solver.solve(system_matrix_pressure, solution_pressure, system_rhs_pressure, preconditioner);
    constraints_pressure.distribute(solution_pressure);
}