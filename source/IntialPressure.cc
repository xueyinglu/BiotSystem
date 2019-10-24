#include "InitialPressure.h"
using namespace std;

double InitialPressure::value(const dealii::Point<dim> &p,
                                    const unsigned int component) const
{

    double t = this->get_time();

    double return_value = 0;

    double x = p(0);
    double y = p(1);

    if (component == 0)
    {
        return_value = cos(x - y);
        //return_value = exp(-x-y*y);
    }
    else if (component == 1)
    {
        return_value = cos(x - y);
        //return_value = exp(-x-y*y);
    }
    else
    {
        cout << "ERROR " << component << endl;
        exit(0);
    }

    return return_value;
}
