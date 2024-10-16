
#include "Molecule.hpp"
#include <cstdio>
#include <iostream>

int main(int argc, char *argv[])
{
	if (argc == 2) {
		Molecule mol = Molecule(string(argv[1]), 0);

		cout << "Number of atoms: " << mol.natoms << endl;
		cout << "Cartesian Coordinates: " << endl;
		mol.print_geom();

		cout << "Interatomic Distances (bohr): " << endl;
		for (int i = 0; i < mol.natoms; i++) {
			for (int j = 0; j < i; j++) {
				printf("%d %d %8.5f\n", i, j, mol.bond(i, j));
			}
		}

		cout << "Bond angles: " << endl;
		for (int i = 0; i < mol.natoms; i++) {
			for (int j = 0; j < mol.natoms; j++) {
				for (int k = 0; k < mol.natoms; k++) {
					if (i != j && j != k && i != k) {
					if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0) {
						printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k)*(180.0/acos(-1)));
						}
					}
				}
			}
		}
	}

	return 0;
}
