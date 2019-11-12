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
    convergence_table.set_tex_caption("L2_u", "L^2-error \\mathbf{u}");

    ofstream error_table_file("error.tex");
    convergence_table.write_tex(error_table_file);

    a_posterior_indicators_table.set_precision("eta_fs", 8);
    a_posterior_indicators_table.set_precision("eta_alg", 8);
    a_posterior_indicators_table.set_precision("eta_time", 8);
    a_posterior_indicators_table.set_precision("eta_flow", 8);
    a_posterior_indicators_table.set_precision("eta_p_residual", 8);
    a_posterior_indicators_table.set_precision("eta_flux_jump", 8);
    a_posterior_indicators_table.set_precision("eta_face_partial_sigma", 8);
    a_posterior_indicators_table.set_precision("eta_face_sigma", 8);
    a_posterior_indicators_table.set_precision("eta_partial_u", 8);
    a_posterior_indicators_table.set_precision("eta_u", 8);
    a_posterior_indicators_table.set_scientific("eta_fs",true);
    a_posterior_indicators_table.set_scientific("eta_alg",true);
    a_posterior_indicators_table.set_scientific("eta_time", true);
    a_posterior_indicators_table.set_scientific("eta_flow", true);
    a_posterior_indicators_table.set_scientific("eta_p_residual", true);
    a_posterior_indicators_table.set_scientific("eta_flux_jump", true);
    a_posterior_indicators_table.set_scientific("eta_face_partial_sigma", true);
    a_posterior_indicators_table.set_scientific("eta_face_sigma", true);
    a_posterior_indicators_table.set_scientific("eta_partial_u", true);
    a_posterior_indicators_table.set_scientific("eta_u", true);
    a_posterior_indicators_table.set_tex_caption("eta_fs", "\\eta_{\\text{fs}}");
    a_posterior_indicators_table.set_tex_caption("eta_alg", "\\eta_{\\text{alg}}");
    a_posterior_indicators_table.set_tex_caption("eta_time", "\\eta_{\\text{time}}");
    a_posterior_indicators_table.set_tex_caption("eta_flow", "\\eta_{\\text{flow}}");
    a_posterior_indicators_table.set_tex_caption("eta_p_residual", "\\eta_{\\text{p\\_residual}}");
    a_posterior_indicators_table.set_tex_caption("eta_flux_jump", "\\eta_{\\text{flux\\_jump}}");
    a_posterior_indicators_table.set_tex_caption("eta_face_partial_sigma", "\\eta_{\\mathcal{E}_{\\partial \\sigma}}");
    a_posterior_indicators_table.set_tex_caption("eta_face_sigma", "\\eta_{\\mathcal{E}_{\\sigma}}");
    a_posterior_indicators_table.set_tex_caption("eta_partial_u", "\\eta_{\\mathcal{T}_{\\partial u}}");
    a_posterior_indicators_table.set_tex_caption("eta_u", "\\eta_{\\mathcal{T}\\_{ u}}");
    ofstream aposterior_table_file("aposteriori.tex");
    a_posterior_indicators_table.write_tex(aposterior_table_file, false);
}