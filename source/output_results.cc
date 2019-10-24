#include "BiotSystem.h"
#include <iostream>
#include <fstream>
using namespace std;

void BiotSystem::output_results() const{
    ofstream output("initial_disp.vtk");
    vector<string> sol_names;
    sol_names.push_back("x_disp");
    sol_names.push_back("y_disp");
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_displacement);
    data_out.add_data_vector(solution_displacement, sol_names);
    data_out.build_patches();
    data_out.write_vtk(output);
}