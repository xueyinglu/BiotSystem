#ifndef AUX_TOOLS_H_
#define AUX_TOOLS_H_
#include <fstream>
#include <sstream>
using namespace std;
using namespace dealii;

class BitmapFile
{
public:
  BitmapFile(const std::string &name);

  double
  get_value(const double x, const double y) const;

private:
  std::vector<double> image_data;
  double hx, hy;
  int nx, ny;

  double
  get_pixel_value(const int i, const int j) const;
};

// The constructor of this class reads in the data that describes
// the obstacle from the given file name.
BitmapFile::BitmapFile(const std::string &name)
    : image_data(0),
      hx(0),
      hy(0),
      nx(0),
      ny(0)
{
  std::ifstream f(name.c_str());
  AssertThrow(f, ExcMessage(std::string("Can't read from file <") +
                            name + ">!"));

  std::string temp;
  getline(f, temp);
  f >> temp;
  if (temp[0] == '#')
    getline(f, temp);

  f >> nx >> ny;

  AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format."));

  for (int k = 0; k < nx * ny; k++)
  {
    unsigned int val;
    f >> val;
    image_data.push_back(val / 255.0);
  }

  hx = 1.0 / (nx - 1);
  hy = 1.0 / (ny - 1);
}

// The following two functions return the value of a given pixel with
// coordinates $i,j$, which we identify with the values of a function
// defined at positions <code>i*hx, j*hy</code>, and at arbitrary
// coordinates $x,y$ where we do a bilinear interpolation between
// point values returned by the first of the two functions. In the
// second function, for each $x,y$, we first compute the (integer)
// location of the nearest pixel coordinate to the bottom left of
// $x,y$, and then compute the coordinates $\xi,\eta$ within this
// pixel. We truncate both kinds of variables from both below
// and above to avoid problems when evaluating the function outside
// of its defined range as may happen due to roundoff errors.
double
BitmapFile::get_pixel_value(const int i,
                            const int j) const
{
  assert(i >= 0 && i < nx);
  assert(j >= 0 && j < ny);
  return image_data[nx * (ny - 1 - j) + i];
}

double
BitmapFile::get_value(const double x,
                      const double y) const
{
  const int ix = std::min(std::max((int)(x / hx), 0), nx - 2);
  const int iy = std::min(std::max((int)(y / hy), 0), ny - 2);

  const double xi = std::min(std::max((x - ix * hx) / hx, 1.), 0.);
  const double eta = std::min(std::max((y - iy * hy) / hy, 1.), 0.);

  return ((1 - xi) * (1 - eta) * get_pixel_value(ix, iy) +
          xi * (1 - eta) * get_pixel_value(ix + 1, iy) +
          (1 - xi) * eta * get_pixel_value(ix, iy + 1) +
          xi * eta * get_pixel_value(ix + 1, iy + 1));
}

template <int dim>
class BitmapFunction : public Function<dim>
{
public:
  BitmapFunction(const std::string &filename,
                 double x1_, double x2_, double y1_, double y2_, double minvalue_, double maxvalue_)
      : Function<dim>(1),
        f(filename), x1(x1_), x2(x2_), y1(y1_), y2(y2_), minvalue(minvalue_), maxvalue(maxvalue_)
  {
  }

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const
  {
    Assert(dim == 2, ExcNotImplemented());
    double x = (p(0) - x1) / (x2 - x1);
    double y = (p(1) - y1) / (y2 - y1);
    return minvalue + f.get_value(x, y) * (maxvalue - minvalue);
  }

private:
  BitmapFile f;
  double x1, x2, y1, y2;
  double minvalue, maxvalue;
};

// Define some tensors for cleaner notation later.
namespace Tensors
{
template <int dim>
inline Tensor<1, dim>
get_grad_p(unsigned int q,
           std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
{
  Tensor<1, dim> grad_p;
  grad_p[0] = old_solution_grads[q][0][0];
  grad_p[1] = old_solution_grads[q][0][1];
  if (dim == 3)
    grad_p[2] = old_solution_grads[q][0][2];

  return grad_p;
}

template <int dim>
inline double
get_deviator_norm(const Tensor<2, dim> deviator)
{

  return std::sqrt(deviator[0][0] * deviator[0][0] + deviator[0][1] * deviator[0][1] + deviator[1][0] * deviator[1][0] + deviator[1][1] * deviator[1][1]);
  if (dim == 3)
  {
    cout << " Wrong get deviator norm: to be fixed" << endl;
  }
}

template <int dim>
inline Tensor<1, dim>
get_grad_pf(
    unsigned int q,
    const std::vector<std::vector<Tensor<1, dim>>> &old_solution_grads)
{
  Tensor<1, dim> grad_pf;
  grad_pf[0] = old_solution_grads[q][dim][0];
  grad_pf[1] = old_solution_grads[q][dim][1];
  if (dim == 3)
    grad_pf[2] = old_solution_grads[q][dim][2];

  return grad_pf;
}

template <int dim>
inline Tensor<2, dim>
get_grad_u(unsigned int q,
           std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
{
  Tensor<2, dim> grad_u;
  grad_u[0][0] = old_solution_grads[q][0][0];
  grad_u[0][1] = old_solution_grads[q][0][1];

  grad_u[1][0] = old_solution_grads[q][1][0];
  grad_u[1][1] = old_solution_grads[q][1][1];
  if (dim == 3)
  {
    grad_u[0][2] = old_solution_grads[q][0][2];

    grad_u[1][2] = old_solution_grads[q][1][2];

    grad_u[2][0] = old_solution_grads[q][2][0];
    grad_u[2][1] = old_solution_grads[q][2][1];
    grad_u[2][2] = old_solution_grads[q][2][2];
  }

  return grad_u;
}

template <int dim>
inline Tensor<2, dim>
get_Identity()
{
  Tensor<2, dim> identity;
  identity[0][0] = 1.0;
  identity[1][1] = 1.0;
  if (dim == 3)
    identity[2][2] = 1.0;

  return identity;
}

template <int dim>
inline Tensor<1, dim>
get_u(unsigned int q,
      std::vector<Vector<double>> old_solution_values)
{
  Tensor<1, dim> u;
  u[0] = old_solution_values[q](0);
  u[1] = old_solution_values[q](1);
  if (dim == 3)
    u[2] = old_solution_values[q](2);

  return u;
}

template <int dim>
inline Tensor<1, dim>
get_u_LinU(const Tensor<1, dim> phi_i_u)
{
  Tensor<1, dim> tmp;
  tmp[0] = phi_i_u[0];
  tmp[1] = phi_i_u[1];
  if (dim == 3)
    tmp[2] = phi_i_u[2];
  return tmp;
}

template <int dim>
inline double
get_divergence_u(const Tensor<2, dim> grad_u)
{
  double tmp;
  if (dim == 2)
  {
    tmp = grad_u[0][0] + grad_u[1][1];
  }
  else if (dim == 3)
  {
    tmp = grad_u[0][0] + grad_u[1][1] + grad_u[2][2];
  }

  return tmp;
}

} // namespace Tensors

// Computing the radius later
template <int dim>
class Function_X : public Function<dim>
{
public:
  Function_X();

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
class Function_Y : public Function<dim>
{
public:
  Function_Y();

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
Function_X<dim>::Function_X() : Function<dim>(1)
{
}

template <int dim>
Function_Y<dim>::Function_Y() : Function<dim>(1)
{
}

// Computing the radius later
template <int dim>
double Function_X<dim>::value(const Point<dim> &p,
                              const unsigned int /*component=0*/) const
{

  double x = p(0);
  double y = p(1);

  double return_value = 0.;

  return x;
}

template <int dim>
double Function_Y<dim>::value(const Point<dim> &p,
                              const unsigned int /*component=0*/) const
{

  double x = p(0);
  double y = p(1);

  double return_value = 0.;

  return y;
}
#endif