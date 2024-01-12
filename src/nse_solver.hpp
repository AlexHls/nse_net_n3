#ifndef NSE_SOLVER_HPP
#define NSE_SOLVER_HPP

#define NSE_MAXITER 500
#define NSE_DIFFVAR 1e-12	/* variation for numerical derivatives */

void run_nse_solver(double temp, double rho, double ye, double *x, struct network_data *nd, struct network_workspace *nw);
void calcprefac(double temp, double rho, struct network_data *nd, struct network_workspace *nw);
void network_part(double temp, const struct network_data *nd, struct network_workspace *nw);

static void guessmu(double ye, double *mu_n, double *mu_p);
static void calcyi(double temp, double rho, double ye, double mu_n, double mu_p,
                   double *yi, double *x, struct network_data *nd,
                   struct network_workspace *nw);

#ifdef NETWORK_SCREENING
static double screening_factor(double Gamma_e, double mu_p_coul, int i, struct network_data *nd, struct network_workspace *nw);
static void screening_params(double rho, double ye, double temp, double *gamma_e, double *mu_p_coul);
#endif /* NETWORK_SCREENING */

#endif // NSE_SOLVER_HPP
