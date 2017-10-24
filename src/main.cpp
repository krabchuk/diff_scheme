#include "cpp_headers.h"
#include "diagonal_matrix.h"
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
  for (int i = 0; i < 5; i++)
    solver.print_residual (0);

  return 0;
}

