#include "DisplacementSolution.h"
using namespace std;
inline Tensor<1, dim> DisplacementSolution::value(const Point<dim> &p) const
{
  Tensor<1, dim> return_tensor;
  double k = 0.05;
  double alpha = 0.75;
  double inv_M = 3. / 28;
  double PI = atan(1) * 4;
  double tx = PI *p(0);
  double ty = PI *p(1);
  double A = 2 * PI * PI * k / (alpha + inv_M);
  return_tensor[0] = -exp(-A * t) / (2 * PI) * cos(tx) * sin(ty);
  return_tensor[1] = -exp(-A * t) / (2 * PI) * sin(tx) * cos(ty);
  return return_tensor;
}

void DisplacementSolution::value_list(const vector<Point<dim>> &points,
                                      vector<Tensor<1, dim>> &value_list) const
{
  Assert(value_list.size() == points.size(),
         ExcDimensionMismatch(value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p = 0; p < n_points; ++p)
    value_list[p] = value(points[p]);
}