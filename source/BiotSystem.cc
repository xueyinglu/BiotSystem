#include "BiotSystem.h"
BiotSystem::BiotSystem() : fe_pressure(1),
                           dof_handler_pressure(triangulation),
                           fe_displacement(FE_Q<dim>(1), dim),
                           dof_handler_displacement(triangulation),
                           permeability(0.05),
                           lambda(0.5),
                           mu(0.125)
{
}
BiotSystem::BiotSystem(int _num_global_refinement, double _del_t, double _T) : fe_pressure(1),
                           dof_handler_pressure(triangulation),
                           fe_displacement(FE_Q<dim>(1), dim),
                           dof_handler_displacement(triangulation),
                           permeability(0.05),
                           lambda(0.5),
                           mu(0.125)
{
    num_global_refinement = _num_global_refinement;
    del_t = _del_t;
    T = _T;
    h = 1./std::pow(2, num_global_refinement);
}