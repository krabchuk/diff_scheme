#include "diagonal_matrix.h"

diagonal_matrix::diagonal_matrix()
{

}

diagonal_matrix::diagonal_matrix (int n_arg)
{
  n = n_arg;

  left = new double [n];
  mid = new double [n];
  right = new double [n];

  rhs = new double [n];
  solution = new double [n];
}

diagonal_matrix::~diagonal_matrix ()
{
  n = 0;

  delete left;
  delete mid;
  delete right;

  delete rhs;
  delete solution;

  left = 0;
  mid = 0;
  right = 0;

  rhs = 0;
  solution = 0;
}
