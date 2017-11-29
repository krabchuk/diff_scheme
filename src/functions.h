#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <math.h>
double rho (double t, double x)
{
  return exp (t) * (x + 1);
//  return t + 1;
//  return t * x;
//  return exp (t);
  return exp(t) * (cos(M_PI * x / 10) + 1.5);
}
double velocity (double t, double x)
{
  return x * (x - 1);
//  return x * (x - 10);
//  return x * (x - 10) * t;
//  return sin (M_PI * x / 5);
  return cos (2 * M_PI * t) * sin (M_PI * x * x / 10 / 10);
}

double h_init (double x)
{
  (void) x;
  return rho(0, x);
}

double v_init (double x)
{
  (void) x;
  return velocity (0, x);
}

// f1 = d(rho)/dt + d(rho * u) / dx
double f1 (double x, double t)
{
//  return  0;
  return exp (t) *(x + 1) + exp (t) * (3 * x * x - 1);
//  return 1 + (t + 1) * (2 * x - 10);
//  return x + t * t * (3 * x * x - 20 * x);
//  return exp (t) * (1 + M_PI / 5 * cos (M_PI * x / 5));
  return rho (x, t) - velocity(x, t) * M_PI / 10 * exp(t) * sin(M_PI * x / 10) + rho(x, t) * cos(2 * M_PI * t) * cos (M_PI * x * x / 100) * M_PI / 50 * x;
}
//f2 = d(rho * u) / dt + d(rho * u * u) / dx + dp/dx - mu d(du/dx)/dx
//p(x) = rho(x) ^ gamma
double f2 (double x, double t, double mu, double gamma)
{
//  return 0;
  return exp (t) * x * (x * x - 1) + exp (t) * (2 * x * (x * x - 1) * (x - 1) + 2 * x * x * x * (x - 1) + x * x * (x * x - 1)) + gamma * exp (t) * pow (rho (x, t), gamma - 1) - mu;
//  return x * (x - 10) +
//  return 2 * t * x * x * (x - 10) + t * t * t * (3 * x * x * (x - 10) * (x - 10) + 2 * x * x * x * (x - 10)) + t - mu * t;
//  return rho (x, t) * velocity (x, t) + rho (x, t) * (2 * sin (M_PI * x / 5) * cos (M_PI * x / 5) * M_PI / 5) + mu * M_PI * M_PI / 5 / 5 * sin (M_PI * x / 5);
  return velocity (x, t) * rho(x, t) - rho(x, t) * 2 * M_PI * sin (2 * M_PI * t) * sin (M_PI * x * x / 100)
      - M_PI / 10 * exp(t) * sin(M_PI * x / 10) * velocity(x, t) * velocity(x, t) + rho(x, t) * 2 * velocity(x, t) * cos(2 * M_PI * t) * cos(M_PI * x * x / 100) * 2 * M_PI * x / 100
      + gamma * pow (rho(x, t), gamma - 1) * (- exp(t) * sin(M_PI * x / 10) * M_PI / 10)
      - mu * cos (2 * M_PI * t) * (-sin(M_PI * x * x / 100) * M_PI * M_PI * x * x / 50 / 50 + cos(M_PI * x * x / 100) * 2 * M_PI / 100);
}

#endif // FUNCTIONS_H
