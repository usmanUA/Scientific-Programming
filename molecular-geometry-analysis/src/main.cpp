
#include "Molecule.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


void	print_bond_angles( int natoms, Molecule& mol)
{
	cout << "\nBond angles: " << endl;
	for (int i = 0; i < natoms; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0) {
					printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k)*(180.0/acos(-1)));
				}
			}
		}
	}
}

void	print_bond_lengths( int natoms, Molecule& mol )
{
	cout << "\nInteratomic Distances (bohr): " << endl;
	for (int i = 0; i < natoms; i++) {
		for (int j = 0; j < i; j++) {
			printf("%d %d %8.5f\n", i, j, mol.bond(i, j));
		}
	}
}

void	print_oop_angles( int natoms, Molecule& mol )
{
	cout << "\nOut of plane Angles: " << endl;
	for (int i = 0; i < natoms; i++) {
		for (int k = 0; k < natoms; k++) {
			for (int j = 0; j < natoms; j++) {
				for (int l = 0; l < j; l++) {
					if (i!=j && i!=k && i!=l && j!=k && k!=l && mol.bond(i,k)<4.0 && mol.bond(k,j)<4.0 && mol.bond(k,l)<4.0) {
						printf("%2d-%2d-%2d-%2d %10.6f\n",i,j,k,l,mol.oop_angle(i,j,k,l)*(180.0/acos(-1.0)));
						}
				}
			}
		}
	}
}

void	print_dihedral_angles( int natoms, Molecule& mol )
{
	cout << "\nDihedral Angles: " << endl;
	for (int i = 0; i < natoms; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				for (int l = 0; l < k; l++) {
					if (mol.bond(i,j)<4.0 && mol.bond(j,k)<4.0 && mol.bond(k,l)<4.0) {
						printf("%2d-%2d-%2d-%2d %10.6f\n",i,j,k,l,mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
						}
				}
			}
		}
	}
}

void	center_mass( int natoms, Molecule& mol )
{
	cout << "\nMolecular center of mass:\t";
	double	x_cm = 0.0;
	double	y_cm = 0.0;
	double	z_cm = 0.0;
	double	M = 0.0;
	double	m_i = 0.0;

	for (int i = 0; i < natoms; i++) {
		m_i = atomic_mass[mol.zvals[i]];
		x_cm += m_i * mol.geom[i][0];
		y_cm += m_i * mol.geom[i][1];
		z_cm += m_i * mol.geom[i][2];
		M += m_i;
	}
	x_cm /= M;
	y_cm /= M;
	z_cm /= M;
	mol.translate(-x_cm, -y_cm, -z_cm);
}
void	principal_inertia( int natoms, Molecule& mol )
{
	double	m_i = 0.0;
	Matrix	I(3,3);

	for (int i = 0; i < natoms; i++) {
		m_i = atomic_mass[mol.zvals[i]];
		I(0,0) += m_i * (pow(mol.geom[i][1],2) + pow(mol.geom[i][2],2));
		I(1,1) += m_i * (pow(mol.geom[i][0],2) + pow(mol.geom[i][2],2));
		I(2,2) += m_i * (pow(mol.geom[i][0],2) + pow(mol.geom[i][1],2));
		I(0,1) -= m_i * mol.geom[i][0] * mol.geom[i][1];
		I(0,2) -= m_i * mol.geom[i][0] * mol.geom[i][2];
		I(1,2) -= m_i * mol.geom[i][1] * mol.geom[i][2];
	}
	I(1,0) = I(0,1), I(2,0) = I(0,2), I(2,1) = I(1,2);
	cout << "\nMoment of Inertia Tensor: (amu * bohr^2)" << endl;
	cout << I << endl;

	// NOTE: principal moments
	Eigen::SelfAdjointEigenSolver<Matrix>	solver(I);
	Matrix	eigenvecs = solver.eigenvectors();
	Matrix	eigenvals = solver.eigenvalues();

	cout << "\nPrincipal Moments of Inertia (amu * AA^w):" << endl;
	cout << eigenvals << endl;
	double	conv = pow(0.529177249,2);
	cout << "\nPrincipal Moments of Inertia (amu * AA^w):" << endl;
	cout << eigenvals * conv << endl;

	conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
	cout << "\nPrincipal Moments of Inertia (amu * AA^w):" << endl;
	cout << eigenvals * conv << endl;

	// NOTE: Rotor
	if (natoms == 2) cout << "\nMolecule is diatomic." << endl;
	else if (eigenvals(0) < 1e-4) cout << "\nMolecule is linear." << endl;
	else if (fabs(eigenvals(0)-fabs(eigenvals(1)))<1e-4 && fabs(eigenvals(1)-fabs(eigenvals(2)))<1e-4)
		cout << "\nMolecule is a spherical top." << endl;
	else if (fabs(eigenvals(0)-fabs(eigenvals(1)))<1e-4 && fabs(eigenvals(1)-fabs(eigenvals(2)))>1e-4)
		cout << "\nMolecule is an oblate symmetric top." << endl;
	else if (fabs(eigenvals(0)-fabs(eigenvals(1)))>1e-4 && fabs(eigenvals(1)-fabs(eigenvals(2)))<1e-4)
		cout << "\nMolecule is a prolate symmetric top." << endl;
	else	cout << "\nMolecule is an asymetric top." << endl;

	// NOTE: rotational constants
	double pi = acos(-1.0);
	double _conv = 6.6260755E-34/(8.0 * pi * pi);
	_conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	conv = _conv * 1e-6;
	cout << "\nRotational constants (MHz):\n";
	cout << "\tA = " << conv/eigenvals(0) << "\t B = " << conv/eigenvals(1) << "\t C = " << conv/eigenvals(2) << endl;
	
	conv = _conv / 2.99792458E10;
	cout << "\nRotational constants (cm-1):\n";
	cout << "\tA = " << conv/eigenvals(0) << "\t B = " << conv/eigenvals(1) << "\t C = " << conv/eigenvals(2) << endl;
}

void	rotational_constants( int natoms, Molecule& mol )
{
}

int main(int argc, char *argv[])
{
	if (argc == 2) {
		Molecule mol = Molecule(string(argv[1]), 0);

		cout << "Number of atoms: " << mol.natoms << endl;
		cout << "Cartesian Coordinates: " << endl;
		mol.print_geom();
		print_bond_lengths(mol.natoms, mol);
		print_bond_angles(mol.natoms, mol);
		print_oop_angles(mol.natoms, mol);
		print_dihedral_angles(mol.natoms, mol);
		center_mass( mol.natoms, mol );
		principal_inertia( mol.natoms, mol );
	}
	return 0;
}
