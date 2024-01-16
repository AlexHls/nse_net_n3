#include <iostream>
#include <stdio.h>
#include <string>

#include "input.hpp"
#include "main.hpp"
#include "nse_solver.hpp"
#include "output.hpp"

int main(int argc, char *argv[]) {
  if (argc != 6) {
    printf("Usage: %s <temp> <rho> <ye> <outfile> <species>\n", argv[0]);
    return 1;
  }

  double temp = std::stod(argv[1]);
  double rho = std::stod(argv[2]);
  double ye = std::stod(argv[3]);
  std::string outfile = argv[4];
  std::string species = argv[5];

  network_data nd;
  network_workspace nw;

  load_data(&nd, species);
  nw.gg = (network_var *)malloc(nd.nuc_count * sizeof(network_var));
  nw.prefac = (double *)malloc(nd.nuc_count * sizeof(double));

  double x[nd.nuc_count];

  run_nse_solver(temp, rho, ye, x, &nd, &nw);

  write_output(outfile, x, &nd);

  return 0;
}
