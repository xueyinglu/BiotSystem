#include "BiotSystem.h"
void BiotSystem::run_fixed_stress(){
    make_grid();
    setup_system();
    assemble_system_displacement();
    solve_displacement();
    output_results();
}