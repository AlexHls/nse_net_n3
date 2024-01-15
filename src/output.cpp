#include "output.hpp"
#include "main.hpp"
#include <iostream>
#include <stdio.h>
#include <string>

void write_output(std::string outfile, double *x, network_data *nd) {
  FILE *fd;
  printf("Writing output...\n");

  if (!(fd = fopen(outfile.c_str(), "w"))) {
    printf("Error opening output file\n");
    exit(1);
  }

  fprintf(fd, "Nucleide, A, Z, Abundance\n");

  for (int i = 0; i < nd->nuc_count; i++) {
    fprintf(fd, "%s, %d, %d, %e\n", nd->nucdata[i].name, nd->nucdata[i].na,
            nd->nucdata[i].nz, x[i]);
  }

  fclose(fd);
}
