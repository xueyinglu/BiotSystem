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