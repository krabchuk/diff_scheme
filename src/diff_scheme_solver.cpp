#include "diff_scheme_solver.h"
#include "cpp_headers.h"

diff_scheme_solver::diff_scheme_solver ()
{
  min_x = 0;
  max_x = 10;

  min_t = 0;
  max_t = 1;

  n = 1;
  m = 1;

  tau = 0;
  h = 0;
  gama = 0;
}

diff_scheme_solver::~diff_scheme_solver ()
{

}

double diff_scheme_solver::ro_0 (double x)
{
  if (x < min_x || x > max_x)
    {
      printf ("ALARM! x = %e\n", x);
      return 0;
    }

  return cos (x * M_PI / 10) + 1.5;
}

double diff_scheme_solver::u_0 (double x)
{
  if (x < min_x || x > max_x)
    {
      printf ("ALARM! x = %e\n", x);
      return 0;
    }

  if (fabs (x - min_x) < MIN_FOR_DIVISION || fabs (x - max_x) < MIN_FOR_DIVISION)
    return 0;

  return sin (x * x * M_PI / 100);
}

double diff_scheme_solver::ro (double t, double x)
{
  if (x < min_x || x > max_x)
    printf ("ALARM! x = %e\n", x);
  if (t < min_t || t > max_t)
    printf ("ALARM! t = %e\n", t);
  return exp (t) * (cos (x * M_PI / 10) + 1.5);
}

double diff_scheme_solver::u (double t, double x)
{
  if (x < min_x || x > max_x)
    printf ("ALARM! x = %e\n", x);
  if (t < min_t || t > max_t)
    printf ("ALARM! t = %e\n", t);
  return cos (2 * M_PI * t) * sin (x * x * M_PI / 100);
}

void diff_scheme_solver::init (int n_arg, int m_arg)
{
  n = n_arg;
  m = m_arg;

  h = (max_x - min_x) / m;
  tau = (max_t - min_t) / n;
  gama = h / tau;

  H_current.resize (m);
  H_next.resize (m);


  V_current.resize (m + 1);
  V_next.resize (m + 1);

  left.resize (m);
  mid.resize (m);
  right.resize (m);
  rhs.resize (m);

  // init first layer
  for (int i = 0; i < m; ++i)
    {
      H_current[i] = ro_0 (min_x + h / 2 + i * h);
      V_current[i] = u_0 (min_x + i * h);
    }
  V_current[m] = 0;
}

void diff_scheme_solver::update_layer ()
{
  for (int i = 0; i < m; ++i)
    {
      H_current[i] = H_next[i];
      V_current[i] = V_next[i];
    }
  V_current[m] = V_next[m];
}

void diff_scheme_solver::build_first_system ()
{
  for (int i = 0; i < m; ++i)
    {
      left[i] = -0.5 * (V_current[i] + fabs (V_current[i]));
      mid[i] = gama + 0.5 * (V_current[i + 1] + fabs (V_current[i + 1]) - V_current[i] + fabs (V_current[i]));
      right[i] = 0.5 * (V_current[i + 1] - fabs (V_current[i + 1]));

      rhs[i] = gama * H_current[i];
    }
}

void diff_scheme_solver::print_residual (int iter)
{
  double t = iter * tau;

  double residual_H = 0;
  double residual_V = 0;

  for (int i = 0; i < m; ++i)
    {
      residual_H += (H_current[i] - ro (t, min_x + h / 2 + i * h)) *
                    (H_current[i] - ro (t, min_x + h / 2 + i * h));
      residual_V += (V_current[i] - u (t, min_x + i * h)) *
                    (V_current[i] - u (t, min_x + i * h));
    }
  residual_H = sqrt (residual_H);
  residual_V = sqrt (residual_V);

  printf ("Iter %d: L2 norm residual H = %e; L2 norm residual V = %e;\n"
          "---------------------------------------------------------------\n",
          iter, residual_H, residual_V);
}

void diff_scheme_solver::solve_system ()
{
  // Solution appear in rhs

  // Forward
  for (int i = 0; i < m - 1; ++i)
    {
      right[i] /= mid[i];
      rhs[i] /= mid[i];

      mid[i + 1] -= right[i] * left[i + 1];
      rhs[i + 1] -= rhs[i] * left[i + 1];
    }

  // Back
  rhs[m - 1] /= mid[m - 1];
  for (int i = m - 2; i >= 0; --i)
    {
      rhs[i] -= rhs[i + 1] * right[i];
    }
}

void diff_scheme_solver::solve_first_system ()
{
  solve_system ();

  for (int i = 0; i < m; ++i)
    {
      H_next[i] = rhs[i];
    }
}

void diff_scheme_solver::print_H ()
{
  for (int i = 0; i < m; ++i)
    {
      printf ("H[%d] = %e\n", i, H_current[i]);
    }
  printf ("----------\n");
}

void diff_scheme_solver::print_V()
{
  for (int i = 0; i < m + 1; ++i)
    {
      printf ("V[%d] = %e\n", i, V_current[i]);
    }
  printf ("----------\n");
}















