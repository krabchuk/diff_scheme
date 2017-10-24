#ifndef DIFF_SCHEME_SOLVER_H
#define DIFF_SCHEME_SOLVER_H

#include "cpp_headers.h"
class diagonal_matrix;

class diff_scheme_solver
{
private:
  double min_x;
  double max_x;

  double min_t;
  double max_t;

  int n;
  int m;

  double tau;
  double h;
  double gama;

  std::vector <double> H_current;
  std::vector <double> H_next;

  std::vector <double> V_current;
  std::vector <double> V_next;

  std::vector <double> left;
  std::vector <double> mid;
  std::vector <double> right;
  std::vector <double> rhs;
public:
  diff_scheme_solver();
  ~diff_scheme_solver ();

  double ro_0 (double x);
  double u_0 (double x);

  double ro (double t, double x);
  double u (double t, double x);

  void init (int n_arg, int m_arg);
  void build_first_system ();

  void update_layer ();

  void print_residual (int iter);
};

#endif // DIFF_SCHEME_SOLVER_H
