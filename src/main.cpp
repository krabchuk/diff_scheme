#include "cpp_headers.h"
#include "diff_scheme_solver.h"

int main(int argc, char *argv [])
{
  if (argc < 2)
    {
      printf ("Usage: %s n m mu\n", argv[0]);
      return 0;
    }

  int n = 0, m = 0;
  double mu = 0;
  if (sscanf (argv[1], "%d", &n) != 1)
    {
      printf ("Usage: %s n m mu\n", argv[0]);
      return 0;
    }
  if (sscanf (argv[2], "%d", &m) != 1)
    {
      printf ("Usage: %s n m mu\n", argv[0]);
      return 0;
    }
  if (sscanf (argv[3], "%lf", &mu) != 1)
    {
      printf ("Usage: %s n m mu\n", argv[0]);
      return 0;
    }

  diff_scheme_solver solver;

  solver.init (n, m, mu);

  solver.fill_debug (0);
  solver.print_residual (0);

  for (int i = 1; i < n; i++)
    {
      solver.calculate_step (i);
//      solver.print_H (i);
//      solver.print_V (i);
    }

  return 0;
}

