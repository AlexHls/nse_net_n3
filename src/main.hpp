#ifndef MAIN_HPP
#define MAIN_HPP

typedef struct {
  double v; /* the values itself */
  double drho, dT, dYe, dYz; /* the derivatives */
  double dabar, dzbar; /* additional derivatives */
} network_var;

struct network_data {
    int nuc_count;
    struct network_nucdata *nucdata;
};

struct network_nucdata {
    int na, nz, nn; /* atomic number, proton number, neutron number */
    char name[6]; /* full name e.g. he4 */
    double part[24]; /* partition function */
    double exm; /* mass excess */
    double m; /* atomic mass */
    double q; /* binding energy */
    double spin; /* spin */
    double nrates, nweakrates; /* nubmer of rates that influence the species */
};

struct network_workspace {
    double *prefac;
    network_var *gg;
};

#endif // MAIN_HPP
