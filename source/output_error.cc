#include "BiotSystem.h"
using namespace std;
#include <iostream>
#include <fstream>
void BiotSystem::output_error() {
    convergence_table.set_precision("L2_p", 4);
    convergence_table.set_precision("L2_u", 4);
    convergence_table.set_scientific("L2_p", true);
    convergence_table.set_scientific("L2_u", true);
    convergence_table.set_tex_caption("L2_p", "L^2-error p");
    convergence_table.set_tex_caption("L2_u", "L^2-error \mathbf{u}");

    ofstream error_table_file("error.tex");
    convergence_table.write_tex(error_table_file);

    a_posterior_indicators_table.set_precision("eta_fs", 8);
    a_posterior_indicators_table.set_precision("eta_alg", 8);
    a_posterior_indicators_table.set_precision("eta_time", 8);
    a_posterior_indicators_table.set_precision("eta_flow", 8);
    a_posterior_indicators_table.set_precision("flux_jump", 8);
    a_posterior_indicators_table.set_scientific("eta_fs",true);
    a_posterior_indicators_table.set_scientific("eta_alg",true);
    a_posterior_indicators_table.set_scientific("eta_time", true);
    a_posterior_indicators_table.set_scientific("eta_flow", true);
    a_posterior_indicators_table.set_scientific("flux_jump", true);
    ofstream aposterior_table_file("aposteriori.tex");
    a_posterior_indicators_table.write_tex(aposterior_table_file);
}