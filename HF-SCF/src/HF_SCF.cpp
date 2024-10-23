

#include "HF_SCF.hpp"
#include <cstdlib>

HFSCF::HFSCF( const char* enuc, string s, string t, string v, string eriout )
{
	S.resize(7,7),T.resize(7,7),V.resize(7,7),ee.resize(7,7),ERI.resize(4000), S_12.resize(7,7);
	FILE*	file = fopen(enuc, "r");
	fscanf(file, "%lf", &this->E_nuc);
	cout << "\033[31mNuclear repulsion energy: \033[0m" << E_nuc << endl;
	cout << "\033[31mOverlap Matrix S: \033[0m" << endl;
	parse_one_e_integrals( s, S );	
	cout << "\033[31mK.E Integrals: \033[0m" << endl;
	parse_one_e_integrals( t, T );	
	cout << "\033[31mP.E Integrals: \033[0m" << endl;
	parse_one_e_integrals( v, V );	
	cout << "\033[31mCore Hamiltonian: \033[0m" << endl;
	compute_core_hamiltonian();
	parse_eri( eriout );
	cout << "\033[31mS_12 : \033[0m" << endl;
	ortogon_S();
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
	cout << mtr << endl << endl;
}

void	HFSCF::compute_core_hamiltonian( void )
{
	H_core = T + V;
	cout << H_core << endl << endl;
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
		ij = i>j ? i*(i+1)/2+j : j*(j+1)/2+i;
		kl = k>l ? k*(k+1)/2+l : k*(k+1)/2+l;
		ijkl = ij>kl ? ij*(ij+1)/2+kl : kl*(kl+1)/2+ij;
		// if (ijkl > 227) {
		// 	cout << "ijkl: " << ijkl << endl;
		// 	exit(1);
		// }
		input >> ERI(ijkl-1);
	}
	// cout << "HERE\n";
	// cout << ERI << endl;
}

void	HFSCF::ortogon_S( void )
{
	Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
	Matrix	eigenvecs = solver.eigenvectors();
	Matrix	eigenvals = solver.eigenvalues();

	Matrix eigen_diags = eigenvals.asDiagonal();
	eigen_diags.resize(7,7);
	//S_12 = eigenvecs * eigenvals * eigenvecs.transpose();
	for (int i = 0; i < eigenvals.size(); i++) {
		eigen_diags(i,i) = 1.0/sqrt(eigenvals(i));
	}
	S_12 = eigenvecs * eigen_diags * eigenvecs.transpose();
	cout <<  S_12 << endl;
}
