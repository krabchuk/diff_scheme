#include "cpp_headers.h"
#include "diff_scheme_solver.h"

int main(int argc, char *argv [])
{
  if (argc < 2)
    {
      printf ("Usage: %s n m", argv[0]);
      return 0;
    }

  int n = 0, m = 0;
  if (sscanf (argv[1], "%d", &n) != 1)
    {
      printf ("Usage: %s n m", argv[0]);
      return 0;
    }
  if (sscanf (argv[2], "%d", &m) != 1)
    {
      printf ("Usage: %s n m", argv[0]);
      return 0;
    }

  diff_scheme_solver solver;

  solver.init (n, m);
  solver.print_H ();
  solver.print_V ();
  solver.build_first_system ();
  solver.solve_first_system ();
  solver.update_layer ();
  solver.print_H ();
  solver.print_residual (1);

  return 0;
}

