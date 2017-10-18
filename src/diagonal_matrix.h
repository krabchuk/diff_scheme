#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "cpp_headers.h"


class diagonal_matrix
{
private:
  int n = 0;

  double *left;
  double *mid;
  double *right;

  double *rhs;

  double *solution;

public:
  diagonal_matrix ();
  diagonal_matrix (int n_arg);

  ~diagonal_matrix ();

  int init_matrix (double *left_arg,
               double *mid_arg,
               double *right_arg,
               double *rhs);

  int solve ();



};

#endif // DIAGONAL_MATRIX_H
