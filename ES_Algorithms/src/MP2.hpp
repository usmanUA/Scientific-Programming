
#pragma once

#include "../Eigen/Dense"
#include "HF_SCF.hpp"

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Array<double, Eigen::Dynamic, 1> Array;

class	MP2
{
	public:
		MP2( Array&, Matrix&, Array& );
		~MP2( void );

		void	Noddy_MO_Transformation( void );
		void	Smarter_MO_Transformation( void );
		void	computeE( void );

	// NOTE: Accessor functions
		double&	E_MP2( void );
	private:
		Matrix&	C;
		Array&	E_MO;
		Array&	AO_ERI;
		Array	MO_ERI;
		double	E_mp2;
};
void	zeros(Matrix& X);
void	mmult(Matrix& A, int a, Matrix& B, int b, Matrix& C, int l, int m, int n);
