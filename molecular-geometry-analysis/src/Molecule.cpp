
#include "Molecule.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

Molecule::Molecule(string filename, int q)
{
	this->charge = q;
	ifstream	input(filename);

	assert(input.good());
	input >> this->natoms;
	cout << "Number of atoms: " << this->natoms << endl;
	this->geom.resize(this->natoms, vector<double>(3));
	this->zvals.resize(this->natoms);
	for (int i = 0; i < this->natoms; i++)
	{
		input >> this->zvals[i] >> this->geom[i][0] >> this->geom[i][1] >> this->geom[i][2];
	}
	input.close();
}

Molecule::~Molecule()
{}

void	Molecule::print_geom()
{
	cout << this->natoms << endl;
	for (unsigned int i = 0; i < this->natoms; i++)
	{
		printf("%d %20.12f %20.12f %20.12f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
	}
}

double	Molecule::unit( int coord, int a, int b)
{
	return	-(this->geom[a][coord] - this->geom[b][coord]) / bond(a, b);
}

void	Molecule::rotate( double phi )
{

}

void	Molecule::translate( double x, double y, double z )
{
	for (int i = 0; i < this->natoms; i++)
	{
		this->geom[i][0] += x;
		this->geom[i][1] += y;
		this->geom[i][2] += z;
	}
}

double	Molecule::bond( int a, int b )
{
	return sqrt(pow(this->geom[a][0] - this->geom[b][0],2)+pow(this->geom[a][1] - this->geom[b][1],2)+pow(this->geom[a][2] - this->geom[b][2],2));
}

double	Molecule::angle( int a, int b, int c )
{
	return acos(unit(0, b, a) * unit(0, b, c) + unit(1, b, a) * unit(1, b, c) + unit(2, b, a) * unit(2, b, c));
}

double	Molecule::torsion( int, int, int, int )
{
	return 1.0;
}

double	Molecule::oop_angle( int a, int b, int c, int d )
{
	double	ebcd_x = unit(1, b, c) * unit(2, b, d) - unit(2, b, c) * unit(1, b, d);
	double	ebcd_y = unit(2, b, c) * unit(0, b, d) - unit(0, b, c) * unit(2, b, d);
	double	ebcd_z = unit(0, b, c) * unit(1, b, d) - unit(1, b, c) * unit(0, b, d);

	double	e_xx = ebcd_x * unit(0, b, a);
	double	e_yy = ebcd_y * unit(1, b, a);
	double	e_zz = ebcd_z * unit(2, b, a);

	double theta = (e_xx+e_yy+e_zz)/sin(angle(b, c, d));
	if (theta < 1.0) theta = asin(-1.0);
	else if (theta > 1.0) theta = asin(1.0);
	else theta = asin(theta);
	return theta;
}
