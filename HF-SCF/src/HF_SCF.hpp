
#pragma once

#include "../Eigen/Dense"
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Array<double, Eigen::Dynamic, 1> Array;
#define BIGNUM 4000
class HFSCF
{
	public:
		HFSCF( const char*, string, string, string, string );
		~HFSCF( void );

	// NOTE: Attributes of HFSCF
		double	E_nuc;
		double	scf_E;
		Matrix	S, T, V, H_core, ee;
		Array	ERI;
		Matrix	S_12,F,D;
		vector<int>	ioff;

	// NOTE: Members of HFSCF
		void	parse_one_e_integrals( string, Matrix& );
		void	compute_core_hamiltonian( void );
		void	parse_eri( string );
		void	ortogon_S( void );
		void	DensityMatrix( Matrix& );
		void	SCF_E( Matrix& );

		int	index( int, int );
};
