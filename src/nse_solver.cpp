#include "main.hpp"
#include <cstring>
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <nse_solver.hpp>
#include <stdio.h>

#define set_zero(x) memset(&(x), 0, sizeof(x))
#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define ELECTRON_CHARGE_ESU (4.80320427e-10)

void run_nse_solver(double temp, double rho, double ye, double *x,
                    network_data *nd, network_workspace *nw) {
  printf("Running NSE solver for temp = %f, rho = %f and ye = %f\n", temp, rho,
         ye);

  int iter, i;
  double mu_n, mu_p, dmu_n, dmu_p;
  double y[2], y2[2];
  double jac[4]; /* order: a11, a12, a21, a22 */
  double det, kt;
  double xsum, nn, ne;
#ifdef NETWORK_SCREENING
  double Gamma_e, mu_p_coul;
#endif // NETWORK_SCREENING

  /* preparation */
  printf("Preparing NSE solver\n");
  calcprefac(temp, rho, nd, nw);
  guessmu(ye, &mu_n, &mu_p);

  /* newton raphson */
  printf("Running Newton-Raphson\n");
  iter = 0;
  while (iter < NSE_MAXITER) {
    calcyi(temp, rho, ye, mu_n, mu_p, y, NULL, nd, nw);
    if (fmax(fabs(y[0]), fabs(y[1])) < 1e-12)
      break;

    dmu_n = fabs(mu_n) * NSE_DIFFVAR;
    if (dmu_n == 0.0)
      dmu_n = NSE_DIFFVAR;
    calcyi(temp, rho, ye, mu_n + dmu_n, mu_p, y2, NULL, nd, nw);
    for (i = 0; i < 2; i++)
      jac[i * 2] = (y2[i] - y[i]) / dmu_n;

    dmu_p = fabs(mu_p) * NSE_DIFFVAR;
    if (dmu_p == 0.0)
      dmu_p = NSE_DIFFVAR;
    calcyi(temp, rho, ye, mu_n, mu_p + dmu_p, y2, NULL, nd, nw);
    for (i = 0; i < 2; i++)
      jac[i * 2 + 1] = (y2[i] - y[i]) / dmu_p;

    det = 1.0 / (jac[0] * jac[3] - jac[1] * jac[2]);
    dmu_n = det * (y[0] * jac[3] - y[1] * jac[1]);
    dmu_p = det * (y[1] * jac[0] - y[0] * jac[2]);

    mu_n -= dmu_n;
    mu_p -= dmu_p;
    iter++;
  }

  kt = 1.0 / (GSL_CONST_CGS_BOLTZMANN * temp);

#ifdef NETWORK_SCREENING
  screening_params(rho, ye, temp, &Gamma_e, &mu_p_coul);
#endif /* NETWORK_SCREENING */

  for (i = 0; i < nd->nuc_count; i++) {
    x[i] =
        nw->prefac[i] * exp(kt * (nd->nucdata[i].nz * mu_p +
                                  nd->nucdata[i].nn * mu_n + nd->nucdata[i].q));
#ifdef NETWORK_SCREENING
    x[i] *= screening_factor(Gamma_e, mu_p_coul, i, nd, nw);
#endif /* NETWORK_SCREENING */
  }
  if (iter < NSE_MAXITER) {
    printf("NSE solver converged after %d iterations\n", iter);
    return;
  } else {
    xsum = 0.0;
    ne = 0.0;
    nn = 0.0;
    for (i = 0; i < nd->nuc_count; i++) {
      xsum += x[i];
      ne += nd->nucdata[i].nz * x[i] / nd->nucdata[i].m;
      nn += nd->nucdata[i].na * x[i] / nd->nucdata[i].m;
    }

    printf("NSE not converged for T=%13.6e, rho=%13.6e, ye=%13.6e, "
           "sum_xi=%13.6e, ye=%13.6e\n",
           temp, rho, ye, xsum, ne / nn);
  }
}

void calcprefac(double temp, double rho, network_data *nd,
                network_workspace *nw) {
  int i;
  network_part(temp, nd, nw);

  for (i = 0; i < nd->nuc_count; i++) {
    nw->prefac[i] =
        nd->nucdata[i].m / rho * (2.0 * nd->nucdata[i].spin + 1.0) *
        nw->gg[i].v *
        pow(2.0 * M_PI * nd->nucdata[i].m * GSL_CONST_CGS_BOLTZMANN * temp /
                (GSL_CONST_CGS_PLANCKS_CONSTANT_H *
                 GSL_CONST_CGS_PLANCKS_CONSTANT_H),
            1.5);
  }
}

void network_part(double temp, const struct network_data *nd,
                  struct network_workspace *nw) {
  /* interpolates partition functions, given the temperature */
  /* \TODO implement partition functions for T > 1e10 K (cf. Rauscher paper?) */
  int index, i;
  double tempLeft, tempRight;
  double dlgLeft, dlgRight;
  double grad;

  static const double network_parttemp[24] = {
      1.0e8, 1.5e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8, 6.0e8, 7.0e8,
      8.0e8, 9.0e8, 1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9, 3.5e9,
      4.0e9, 4.5e9, 5.0e9, 6.0e9, 7.0e9, 8.0e9, 9.0e9, 1.0e10};

  index = 0;
  temp = min(max(temp, network_parttemp[0]), network_parttemp[23]);

  while (temp > network_parttemp[index]) {
    index++;
  }
  if (index > 0)
    index--;

  tempLeft = network_parttemp[index];
  tempRight = network_parttemp[index + 1];

  for (i = 0; i < nd->nuc_count; i++) {
    dlgLeft = nd->nucdata[i].part[index];
    dlgRight = nd->nucdata[i].part[index + 1];

    grad = (dlgRight - dlgLeft) / (tempRight - tempLeft);
    set_zero(nw->gg[i]);
    nw->gg[i].v = exp(dlgLeft + (temp - tempLeft) * grad);
    nw->gg[i].dT = nw->gg[i].v * grad;
  }
}

static void guessmu(double ye, double *mu_n, double *mu_p) {
  if (ye > 0.55) {
    *mu_n = -2.7e-5;
    *mu_p = -1.1e-6;
  } else if (ye > 0.53) {
    *mu_n = -0.22e-4;
    *mu_p = -0.06e-4;
  } else if (ye > 0.46) {
    *mu_n = -0.19e-4;
    *mu_p = -0.08e-4;
  } else if (ye > 0.44) {
    *mu_n = -0.145e-4;
    *mu_p = -0.135e-4;
  } else if (ye > 0.42) {
    *mu_n = -0.11e-4;
    *mu_p = -0.18e-4;
  } else {
    *mu_n = -0.08e-4;
    *mu_p = -0.22e-4;
  }
}

static void calcyi(double temp, double rho, double ye, double mu_n, double mu_p,
                   double *yi, double *x, struct network_data *nd,
                   struct network_workspace *nw) {
  int i;
  double kt, xi;
#ifdef NETWORK_SCREENING
  double Gamma_e, mu_p_coul;
#endif /* NETWORK_SCREENING */

  yi[0] = -1.0;
  yi[1] = 0.0;

  kt = 1.0 / (GSL_CONST_CGS_BOLTZMANN * temp);
#ifdef NETWORK_SCREENING
  screening_params(rho, ye, temp, &Gamma_e, &mu_p_coul);
#endif /* NETWORK_SCREENING */
  for (i = 0; i < nd->nuc_count; i++) {
    xi =
        nw->prefac[i] * exp(kt * (nd->nucdata[i].nz * mu_p +
                                  nd->nucdata[i].nn * mu_n + nd->nucdata[i].q));
#ifdef NETWORK_SCREENING
    xi *= screening_factor(Gamma_e, mu_p_coul, i, nd, nw);
#endif /* NETWORK_SCREENING */
    yi[0] += xi;
    yi[1] += GSL_CONST_CGS_UNIFIED_ATOMIC_MASS / nd->nucdata[i].m *
             ((ye - 1) * nd->nucdata[i].nz + ye * nd->nucdata[i].nn) * xi;

    if (x)
      x[i] = xi;
  }
}

#ifdef NETWORK_SCREENING
static double screening_factor(double Gamma_e, double mu_p_coul, int i,
                               struct network_data *nd,
                               struct network_workspace *nw) {
  /* screening as mentioned in Seitenzahl et al. Atomic and Nuclear Data Tables
   * 95 (2009) 96-114 eq. (14) */
  const double A1 = -0.9052, A2 = 0.6322, A3 = -sqrt(3.0) / 2.0 - A1 / sqrt(A2);
  double Gamma_i;

  if (nd->nucdata[i].nz != 0.0) {
    Gamma_i = Gamma_e * pow(nd->nucdata[i].nz, 5.0 / 3.0);
    return exp(
        -(A1 * (sqrt(Gamma_i * (A2 + Gamma_i)) -
                A2 * log(sqrt(Gamma_i / A2) + sqrt(1.0 + Gamma_i / A2))) +
          2.0 * A3 * (sqrt(Gamma_i) - atan(sqrt(Gamma_i)))) +
        nd->nucdata[i].nz * mu_p_coul);
  } else
    return 1.0;
}

static void screening_params(double rho, double ye, double T, double *Gamma_e,
                             double *mu_p_coul) {
  const double A1 = -0.9052, A2 = 0.6322, A3 = -sqrt(3.0) / 2.0 - A1 / sqrt(A2);

  *Gamma_e =
      gsl_pow_2(ELECTRON_CHARGE_ESU) / GSL_CONST_CGS_BOLTZMANN / T *
      pow(4.0 / 3.0 * M_PI * rho * ye / GSL_CONST_CGS_UNIFIED_ATOMIC_MASS,
          1.0 / 3.0);
  *mu_p_coul =
      (A1 * (sqrt(*Gamma_e * (A2 + *Gamma_e)) -
             A2 * log(sqrt(*Gamma_e / A2) + sqrt(1.0 + *Gamma_e / A2))) +
       2.0 * A3 * (sqrt(*Gamma_e) - atan(sqrt(*Gamma_e))));
}
#endif /* NETWORK_SCREENING */
