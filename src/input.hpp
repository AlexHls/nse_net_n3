#ifndef INPUT_HPP
#define INPUT_HPP
#include <string>

void load_data(struct network_data *nd, std::string species);
void read_species(std::string speciesfile, struct network_data *nd);
void read_mass(std::string massfile, struct network_data *nd);
void read_part(std::string partfile, struct network_data *nd);

#endif // INPUT_HPP
