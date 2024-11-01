
#pragma once

#include "../Eigen/Dense"
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Array<double, Eigen::Dynamic, 1> Array;
#define BIGNUM 4000

class SCF
{
	public:
		SCF( const char*, char*, char*, char*, char* );
		~SCF( void );


	// NOTE: Computations
		void	initialize( void );
		void	runSCF( double, double );

	// NOTE: Accessor functions
		Matrix&	AOFock( void );
		Matrix&	DensityMatrix( void );
		Matrix&	CoefficientMatrix( void );
		Array&	MO_Energies( void );
		Array&	AOERI( void );
		double&	E( void );

	private:
	// NOTE: Attributes of SCF
		double	E_nuc;
		double	E_SCF;
		Matrix	S,T,V,H_core,C;
		Array	ERI,E_MOs;
		Matrix	S_12,F,D,AOF;
	// NOTE: Computations
		void	computeD( bool );
		void	computeE( void );
		void	computeAOF( void );
};

void	parse_eri( char *file, Array& ERI );
void	parse_1e_integrals( char* file, Matrix& mtr );
void	zero_out( Matrix& m );
int	index( int i, int j );
void	ioff_matrix( void );
