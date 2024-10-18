
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
	printf("%10.6f\t%10.6f\t%10.6f\n", x, y, z);
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

double	Molecule::torsion( int a, int b, int c, int d )
{
	double	eab_x = unit(0, b, a);
	double	eab_y = unit(1, b, a);
	double	eab_z = unit(2, b, a);

	double	ebc_x = unit(0, c, b);
	double	ebc_y = unit(1, c, b);
	double	ebc_z = unit(2, c, b);

	double	ecd_x = unit(0, d, c);
	double	ecd_y = unit(1, d, c);
	double	ecd_z = unit(2, d, c);

	double	eabc_x = eab_y * ebc_z - eab_z * ebc_y;
	double	eabc_y = eab_z * ebc_x - eab_x * ebc_z;
	double	eabc_z = eab_x * ebc_y - eab_y * ebc_x;

	double	ebcd_x = ebc_y * ecd_z - ebc_z * ecd_y;
	double	ebcd_y = ebc_z * ecd_x - ebc_x * ecd_z;
	double	ebcd_z = ebc_x * ecd_y - ebc_y * ecd_x;

	double	exx = eabc_x * ebcd_x;
	double	eyy = eabc_y * ebcd_y;
	double	ezz = eabc_z * ebcd_z;
	double tau = (exx+eyy+ezz)/(sin(angle(a, b, c)) * sin(angle(b, c, d)));
	if (tau < -1.0) tau = acos(-1.0);
	else if (tau > 1.0) tau = acos(1.0);
	else tau = acos(tau);

	// NOTE: sign of the angle
	double	cross_term_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
	double	cross_term_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
	double	cross_term_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
	double	norm = cross_term_x*cross_term_x + cross_term_y*cross_term_y + cross_term_z*cross_term_z;
	cross_term_x /= norm;
	cross_term_y /= norm;
	cross_term_z /= norm;
	double dot = cross_term_x*unit(0, b, c) + cross_term_y*unit(1, b, c) + cross_term_z*unit(2, b, c);
	double	sign = 1.0;
	if (dot < 0.0) sign = -1.0;
	return tau * sign;
}

double	Molecule::oop_angle( int a, int b, int c, int d )
{
	double	ebcd_x = unit(1, c, b) * unit(2, c, d) - unit(2, c, b) * unit(1, c, d);
	double	ebcd_y = unit(2, c, b) * unit(0, c, d) - unit(0, c, b) * unit(2, c, d);
	double	ebcd_z = unit(0, c, b) * unit(1, c, d) - unit(1, c, b) * unit(0, c, d);

	double	e_xx = ebcd_x * unit(0, c, a);
	double	e_yy = ebcd_y * unit(1, c, a);
	double	e_zz = ebcd_z * unit(2, c, a);

	double theta = (e_xx+e_yy+e_zz)/sin(angle(b, c, d));
	if (theta < -1.0) theta = asin(-1.0);
	else if (theta > 1.0) theta = asin(1.0);
	else theta = asin(theta);
	return theta;
}

