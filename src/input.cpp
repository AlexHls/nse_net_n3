#include "input.hpp"
#include "main.hpp"
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_num.h>
#include <math.h>
#include <stdio.h>
#include <string>

static const double conv =
    1.602177e-12 * 1.0e3 * 6.0221367e23; /* eV2erg * 1.0e3 [keV] * avogadro */

void load_data(network_data *nd, std::string species) {
  printf("Loading data\n");

  std::string massfile = "data/mass.txt";
  std::string partfile = "data/part.txt";

  read_species(species, nd);
  printf("Found %d species\n", nd->nuc_count);

  read_mass(massfile, nd);

  read_part(partfile, nd);
}

void read_species(std::string speciesfile, network_data *nd) {
  FILE *fd;
  char cdummy[100];
  printf("Reading species\n");

  if (!(fd = fopen(speciesfile.c_str(), "r"))) {
    printf("Error opening species file\n");
    exit(1);
  }
  fgets(cdummy, 100, fd);
  sscanf(cdummy, "%u", &nd->nuc_count);

  /* allocate memory for nucdata */
  nd->nucdata =
      (network_nucdata *)malloc(nd->nuc_count * sizeof(network_nucdata));

  for (int i = 0; i < nd->nuc_count; i++) {
    fgets(cdummy, 100, fd);
    sscanf(cdummy, "%5c%d%d", nd->nucdata[i].name, &nd->nucdata[i].na,
           &nd->nucdata[i].nz);
    nd->nucdata[i].name[5] = 0;
    nd->nucdata[i].nn = nd->nucdata[i].na - nd->nucdata[i].nz;
    nd->nucdata[i].nrates = 0;
    nd->nucdata[i].nweakrates = 0;
    nd->nucdata[i].spin = 0.0;
  }
  fclose(fd);
}

void read_mass(std::string massfile, network_data *nd) {
  FILE *fd;
  char cdummy[200];
  int nminz, nn, nz, na;
  double exm, q;

  printf("Reading mass\n");
  if (!(fd = fopen(massfile.c_str(), "r"))) {
    printf("Error opening mass file\n");
    exit(1);
  }

  int *masses;
  masses = (int *)malloc(nd->nuc_count * sizeof(int));
  for (int i = 0; i < nd->nuc_count; i++) {
    masses[i] = 0;
  }

  /* skip 39 lines */
  for (int i = 0; i < 39; i++) {
    fgets(cdummy, 200, fd);
  }

  while (!feof(fd)) {
    if (fgets(cdummy, 200, fd) == NULL) {
      break;
    }
    sscanf(&cdummy[1], "%d%d%d%d\n", &nminz, &nn, &nz, &na);

    for (int i = 0; i < nd->nuc_count; i++) {
      if (nd->nucdata[i].na == na && nd->nucdata[i].nz == nz) {
        sscanf(&cdummy[29], "%lf", &exm);
        nd->nucdata[i].exm = exm;
        sscanf(&cdummy[54], "%lf", &q);
        nd->nucdata[i].q =
            q * 1.602177e-12 * 1.0e3 *
            nd->nucdata[i]
                .na; /* values are in KeV / nucleon, converting to erg */
        masses[i] = 1;
        break;
      }
    }
  }

  fclose(fd);

  int missing = 0;
  for (int i = 0; i < nd->nuc_count; i++) {
    if (masses[i] == 0) {
      printf("Nucleus %s missing in mass file.\n", nd->nucdata[i].name);
      missing += 1;
    }

    nd->nucdata[i].m =
        nd->nucdata[i].na * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS +
        nd->nucdata[i].exm * conv / GSL_CONST_NUM_AVOGADRO /
            (GSL_CONST_CGS_SPEED_OF_LIGHT * GSL_CONST_CGS_SPEED_OF_LIGHT);
  }

  if (missing > 0) {
    exit(1);
  }

  free(masses);
}

void read_part(std::string partfile, network_data *nd) {
  FILE *fd;
  char cdummy[200];
  double spin;
  int nz, na, found;

  printf("Reading partition function\n");

  if (!(fd = fopen(partfile.c_str(), "r"))) {
    printf("Error opening partition file\n");
    exit(1);
  }

  int *spins;
  spins = (int *)malloc(nd->nuc_count * sizeof(int));
  for (int i = 0; i < nd->nuc_count; i++) {
    spins[i] = 0;
  }

  /* default value for the partition function is 1.0 => log(part) = 0.0 */
  for (int i = 0; i < nd->nuc_count; i++) {
    for (int j = 0; j < 24; j++) {
      nd->nucdata[i].part[j] = 0.0;
    }
  }

  /* skip 4 lines */
  for (int i = 0; i < 4; i++) {
    fgets(cdummy, 200, fd);
  }

  while (!feof(fd)) {
    if (fgets(cdummy, 200, fd) == NULL)
      break; /* skip name of the nucleus */
    if (fgets(cdummy, 200, fd) == NULL)
      break;
    sscanf(cdummy, "%d%d%lf\n", &nz, &na, &spin);

    found = 0;
    for (int i = 0; i < nd->nuc_count; i++) {
      if (nd->nucdata[i].na == na && nd->nucdata[i].nz == nz) {
        found = 1;

        /* read partition function data */
        for (int j = 0; j < 3; j++) {
          char *tmp, *next_field;

          fgets(cdummy, 200, fd);
          next_field = cdummy;
          for (int k = 0; k < 8; k++) {
            nd->nucdata[i].part[j * 8 + k] = strtod(next_field, &tmp);
            if (tmp == next_field) {
              printf("Error reading partition function data for element %s\n",
                     nd->nucdata[i].name);
              exit(1);
            }
            next_field = tmp;
            nd->nucdata[i].part[j * 8 + k] =
                log(nd->nucdata[i].part[j * 8 + k]);
          }
        }

        nd->nucdata[i].spin = spin;
        spins[i] = 1;
        break;
      }
    }

    /* skip 3 lines containing partition function data */
    if (!found) {
      for (int i = 0; i < 3; i++) {
        fgets(cdummy, 200, fd);
      }
    }
  }

  fclose(fd);

  for (int i = 0; i < nd->nuc_count; i++) {
    if (spins[i] == 0) {
      printf("There are no partition function and spin data for nucleus %s. "
             "Assuming spin 0 and constant partition function of 1.\n",
             nd->nucdata[i].name);
    }
  }
  free(spins);
}
