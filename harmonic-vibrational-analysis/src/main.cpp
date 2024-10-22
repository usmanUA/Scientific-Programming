
#include "Molecule.hpp"
#include <iostream>
#include <cmath>
#include <cstdio>

int main(int argc, char *argv[])
{
	if (argc == 3) {
		Molecule mol = Molecule(string(argv[1]), 0);

		cout << "Cartesian Coordinates: " << endl;
		mol.print_geom();
		mol.save_hessian(argv[2]);
		mol.mass_weight_hessian();
		mol.diagonalize_mw_hessian();
		cout << endl;
		mol.harmonic_vibrational_frequencies();
	}
	return 0;
}
