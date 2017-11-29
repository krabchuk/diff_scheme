#include "diff_scheme_solver.h"
#include "cpp_headers.h"
#include "functions.h"

diff_scheme_solver::diff_scheme_solver ()
{
  min_x = 0;
  max_x = 1;

  min_t = 0;
  max_t = 1;

  n = 1;
  m = 1;

  tau = 0;
  h = 0;
  gamma = 0.25;

  mu = 0;
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

void diff_scheme_solver::init (int n_arg, int m_arg, double mu_arg)
{
  n = n_arg;
  m = m_arg;

  h = (max_x - min_x) / m;
  tau = (max_t - min_t) / n;
  gamma = 0.25;

  mu = mu_arg;

  H_current.resize (m);
  H_next.resize (m);


  V_current.resize (m + 1);
  V_next.resize (m + 1);

  left.resize (m);
  mid.resize (m);
  right.resize (m);
  rhs.resize (m);

  debug_vector1.resize (m);
  debug_vector2.resize (m + 1);
//  for (int i = 0; i < m; ++i)
//    debug_vector[i] = 0;

  // init first layer
  for (int i = 0; i < m; ++i)
    {
      H_current[i] = h_init (min_x + h / 2 + i * h);
      V_current[i] = v_init (min_x + i * h);
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
      left[i] = -tau / h * 0.5 * (V_current[i] + fabs (V_current[i]));
      mid[i] = 1 + tau / h * 0.5 * (V_current[i + 1] + fabs (V_current[i + 1]) - V_current[i] + fabs (V_current[i]));
      right[i] = tau / h * 0.5 * (V_current[i + 1] - fabs (V_current[i + 1]));

      rhs[i] = H_current[i] + tau * debug_vector1[i];
    }
}

void diff_scheme_solver::build_second_system ()
{
  left[0] = 0;
  mid[0] = 1;
  right[0] = 0;
  rhs[0] = 0;

  for (int i = 1; i < m; ++i)
    {
      double C = 0.5 * (H_next[i] + H_next[i - 1]);

      if (fabs (C) < MIN_FOR_DIVISION)
        {
          left[i] = 0;
          mid[i] = 1;
          right[i] = 0;
          rhs[i] = 0;
        }
      else
        {
          double D = gamma * (pow(fabs (H_next[i]), gamma - 1) - pow(fabs (H_next[i - 1]), gamma - 1)) / (gamma - 1) / h;

          left[i] = - 0.5 * C * (V_current[i] + fabs (V_current[i])) * tau / h - mu * tau / h / h;
          mid[i] = C  + C * fabs (V_current[i]) * tau / h + 2 * mu * tau / h / h;
          right[i] = 0.5 * C * (V_current[i] - fabs (V_current[i])) * tau / h - mu * tau / h / h;
          rhs[i] = C * V_current[i] * tau - C * D * tau + tau * debug_vector2[i];
        }
    }

  left[m] = 0;
  mid[m] = 1;
  right[m] = 0;
  rhs[m] = 0;
}

void diff_scheme_solver::print_residual (int iter)
{
  double t = iter * tau;

  double residual_H = 0;
  double residual_V = 0;
  double max_H = 0;
  double max_V = 0;
  double border_H = 0;
  double border_V = 0;

  for (int i = 0; i < m; ++i)
    {
      residual_H += (H_current[i] - rho (t, min_x + h / 2 + i * h)) *
                    (H_current[i] - rho (t, min_x + h / 2 + i * h));
      residual_V += (V_current[i] - velocity (t, min_x + i * h)) *
                    (V_current[i] - velocity (t, min_x + i * h));

      if (fabs (H_current[i] - rho (t, min_x + h / 2 + i * h)) > max_H)
        max_H = fabs (H_current[i] - rho (t, min_x + h / 2 + i * h));

      if (fabs (V_current[i] - velocity (t, min_x + i * h)) > max_V)
        max_V = fabs (V_current[i] - velocity (t, min_x + i * h));

      if (i == 0 || i == m - 1)
        {
          border_H += residual_H;
          border_V += residual_V;
        }
    }
  residual_H *= h;
  residual_V *= h;
  border_H *= h;
  border_V *= h;


  // border

  printf ("============================Iter %d: ==========================\n"
          "H L2-residual = %e; V L2-residual = %e;\n"
          "H  C-residual = %e; V  C-residual = %e;\n"
          "H  2-residual = %e; V  2-residual = %e;\n"
          "---------------------------------------------------------------\n",
          iter, sqrt (residual_H), sqrt (residual_V),
          max_H, max_V,
          sqrt (sqrt (residual_H) + 0.5 * border_H), sqrt (sqrt (residual_H) + 0.5 * border_H));
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

void diff_scheme_solver::solve_second_system ()
{
  solve_system ();

  for (int i = 0; i < m + 1; ++i)
    {
      V_next[i] = rhs[i];
    }
}

void diff_scheme_solver::print_H (int iter)
{
  for (int i = 0; i < m && i < 10; ++i)
    {
      printf ("H[%d] = %e %e\n", i, H_current[i], rho(iter * tau, min_x + h / 2 + i * h));
    }
  printf ("----------\n");
}

void diff_scheme_solver::print_V(int iter)
{
  for (int i = 0; i < m + 1 && i < 10; ++i)
    {
      printf ("V[%d] = %e %e\n", i, V_current[i], velocity (iter * tau, min_x + i * h));
    }
  printf ("----------\n");
}

void diff_scheme_solver::calculate_step (int iter)
{

  build_first_system ();
  solve_first_system ();

  build_second_system ();
  solve_second_system ();

  update_layer ();
  print_residual (iter);

  fill_debug (iter);
}

void diff_scheme_solver::fill_debug (int iter)
{

  for (int i = 0; i <= m; i++)
    {
      if (i < m)
        debug_vector1[i] = f1 (min_x + i * h + h / 2, tau * iter);
      debug_vector2[i] = f2 (min_x + i * h, tau * iter, mu, gamma);
    }
}













