#ifndef DISPLACEMENT_SOLUTION_H_
#define DISPLACEMENT_SOLUTION_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;

class DisplacementSolution {
    public :
    DisplacementSolution( double _t){t = _t;};
    
    virtual Tensor<1,dim> value(const Point<dim> & point) const;
    virtual void value_list(const vector<Point<dim> > & points, vector<Tensor<1,dim> > & value_list) const ;
    private :
    double t;
};

#endif