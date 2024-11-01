#include "HF_SCF.hpp"

SCF::SCF( const char* enuc, char* s, char* t, char* v, char* eriout )
{
	S.resize(7,7),T.resize(7,7),V.resize(7,7),ERI.resize(BIGNUM), S_12.resize(7,7);
	D.resize(7,7),C.resize(7,7),F.resize(7,7),AOF.resize(7,7);
	FILE*	file = fopen(enuc, "r");
	fscanf(file, "%lf", &E_nuc);
	parse_1e_integrals( s, S );
	parse_1e_integrals( t, T );
	parse_1e_integrals( v, V );
	parse_eri( eriout, ERI );
}

SCF::~SCF( void )
{}

void	SCF::initialize( void )
{
	H_core = T + V;

	Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
	Matrix	eigenvecs = solver.eigenvectors();
	Matrix	eigenvals = solver.eigenvalues();

	Matrix eigen_diags = eigenvals.asDiagonal();
	for (int i = 0; i < eigenvals.size(); i++) {
		eigen_diags(i,i) = 1.0/sqrt(eigenvals(i));
	}
	S_12 = eigenvecs * eigen_diags * eigenvecs.transpose();
}

void	SCF::runSCF( double delta1, double delta2 )
{
	Matrix	prev_D;
	double	prev_E,delta_E,rms;
	int	iter = 0;
	computeD( true );
	computeAOF();
	computeE();
	cout << "HERE\n";
	printf("Iter   \tE(elec)   \tE(tot)   \tDelta(E)   \tRMS(D)\n");
	printf("%d  %12.12f  %12.12f  \n",iter,E_SCF-E_nuc,E_SCF);
	while (true) {
		iter++;
		prev_E = E_SCF;
		prev_D = D;
		computeAOF();
		computeD(false);
		computeE();
		delta_E = E_SCF - prev_E;
		rms = sqrt(((D-prev_D)*(D-prev_D)).sum());
		printf("%d  %12.12f  %12.12f  %12.12f  %12.12f\n",iter,E_SCF-E_nuc,E_SCF,delta_E,rms);
		if (delta_E<delta1 && rms<delta2) {
			break;
		}
	}
}

void	SCF::computeD( bool core )
{
	if (core) {
		F = S_12.transpose() * H_core * S_12;
	}
	else {
		F = S_12.transpose() * AOF * S_12;
	}
	Eigen::SelfAdjointEigenSolver<Matrix> solver(F);
	C = solver.eigenvectors();
	Matrix eigenvals = solver.eigenvalues();
	E_MOs.resize(eigenvals.size());
	E_MOs = eigenvals;
	C = S_12 * C;
	Matrix C_occ = C.leftCols(5);
	D = C_occ * C_occ.transpose();
}

void	SCF::computeE( void )
{
	Matrix	H_F = H_core + AOF;
	E_SCF = (D * H_F).sum() + E_nuc;
}

void	SCF::computeAOF( void )
{
	int ij,kl,ijkl,ik,jl,ikjl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			AOF(i,j) = H_core(i,j);
			for (int k = 0; k < 7; k++) {
				for (int l = 0; l < 7; l++) {
					ij = index(i,j);
					kl = index(k,l);
					ijkl = index(ij,kl);
					ik = index(i, k);
					jl = index(j, l);
					ikjl = index(ik, jl);
					AOF(i,j) += D(k,l) * (2*ERI(ijkl) - ERI(ikjl));
				}
			}
		}
	}
}

Matrix&	SCF::DensityMatrix( void )
{
	return D;
}

Matrix&	SCF::CoefficientMatrix( void )
{
	return C;
}

Matrix&	SCF::AOFock( void )
{
	return AOF;
}

Array&	SCF::MO_Energies( void )
{
	return E_MOs;
}

Array&	SCF::AOERI( void )
{
	return ERI;
}

double&	SCF::E( void )
{
	return E_SCF;
}
