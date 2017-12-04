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
  double gamma;

  double mu;

  std::vector <double> H_current;
  std::vector <double> H_next;

  std::vector <double> V_current;
  std::vector <double> V_next;

  std::vector <double> left;
  std::vector <double> mid;
  std::vector <double> right;
  std::vector <double> rhs;

  std::vector <double> debug_vector1;
  std::vector <double> debug_vector2;
public:
  diff_scheme_solver();
  ~diff_scheme_solver ();

  double ro_0 (double x);
  double u_0 (double x);

  double ro (double t, double x);
  double u (double t, double x);

  void init (int n_arg, int m_arg, double mu_arg);
  void build_first_system ();
  void build_second_system ();

  void update_layer ();

  void print_residual (int iter);

  void solve_first_system ();
  void solve_second_system ();

  void solve_system (int size);

  void print_H (int iter);
  void print_V (int iter);

  void calculate_step (int iter);

  void fill_debug (int iter);
};

#endif // DIFF_SCHEME_SOLVER_H
