

#include "HF_SCF.hpp"
#include <cstdlib>
#include <vector>

HFSCF::HFSCF( const char* enuc, string s, string t, string v, string eriout )
{
	ioff.resize(BIGNUM);
	ioff[0] = 0;
	for (int i = 1; i < BIGNUM; i++) {
		ioff[i] = ioff[i-1] + i;
	}
	S.resize(7,7),T.resize(7,7),V.resize(7,7),ee.resize(7,7),ERI.resize(BIGNUM), S_12.resize(7,7);
	D.resize(7,7),F_bar.resize(7,7),F.resize(7,7);
	FILE*	file = fopen(enuc, "r");
	fscanf(file, "%lf", &E_nuc);
	parse_one_e_integrals( s, S );	
	parse_one_e_integrals( t, T );	
	parse_one_e_integrals( v, V );	
	parse_eri( eriout );
}

HFSCF::~HFSCF( void )
{}

void	HFSCF::parse_one_e_integrals( string file, Matrix& mtr )
{
	ifstream	input(file);
	int x,y;
	if (!input.is_open()) {
		cout << "Unable to open: " << file << endl;
		exit(1);
	}
	while (input.peek() != EOF) {
		input >> x >> y >> mtr(x-1,y-1);
		mtr(y-1,x-1) = mtr(x-1,y-1);
	}
}

void	HFSCF::compute_core_hamiltonian( void )
{
	H_core = T + V;
}

void	HFSCF::parse_eri( string file )
{
	ifstream	input(file);
	int i,j,k,l,ij,kl,ijkl;
	if (!input.is_open()) {
		cout << "Unable to open: " << file << endl;
		exit(1);
	}
	while (input.peek() != EOF) {
		input >> i >> j >> k >> l;
		ij = index(i, j);
		kl = index(k, l);
		ijkl = index(ij, kl);
		input >> ERI(ijkl);
	}
}

void	HFSCF::ortogon_S( void )
{
	Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
	Matrix	eigenvecs = solver.eigenvectors();
	Matrix	eigenvals = solver.eigenvalues();

	Matrix eigen_diags = eigenvals.asDiagonal();
	eigen_diags.resize(7,7);
	for (int i = 0; i < eigenvals.size(); i++) {
		eigen_diags(i,i) = 1.0/sqrt(eigenvals(i));
	}
	S_12 = eigenvecs * eigen_diags * eigenvecs.transpose();
	for (int i = 0; i < S_12.size()/7; i++) {
		for (int j = 0; j < S_12.size()/7; j++) {
			if (abs(S_12(i,j)) < 1e-10){
				S_12(i,j) = 0.00000000;
			}
		}

	}
}

void	HFSCF::DensityMatrix( Matrix& F_matrix )
{
	F_bar = S_12.transpose() * F_matrix * S_12;
	for (int i = 0; i < F_bar.size()/7; i++) {
		for (int j = 0; j < F_bar.size()/7; j++) {
			if (abs(F_bar(i,j)) < 1e-10){
				F_bar(i,j) = 0.00000000;
			}
		}

	}
	Eigen::SelfAdjointEigenSolver<Matrix> solver(F_bar);
	Matrix eigenvecs = solver.eigenvectors();
	Matrix eigenvals = solver.eigenvalues();
	eigenvecs = S_12 * eigenvecs;
	for (int i = 0; i < eigenvecs.size()/7; i++) {
		for (int j = 0; j < eigenvecs.size()/7; j++) {
			if (abs(eigenvecs(i,j)) < 1e-10){
				eigenvecs(i,j) = 0.00000000;
			}
		}

	}
	D = eigenvecs * eigenvecs.transpose();
}

void	HFSCF::SCF_E( Matrix& F_uv )
{
	Matrix	H_F = H_core + F_uv;
	scf_E = (D * H_F.transpose()).sum() + E_nuc;
}

int	HFSCF::index( int i, int j ) 
{
	int ij = i>j ? ioff[i]+j : ioff[j]+i;
	return ij;
}
