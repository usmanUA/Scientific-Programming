
#include "HF_SCF.hpp"

void	new_fock( HFSCF& scf, Matrix& F )
{
///	cout << "DensityMatrix: \n" << scf.D << endl;
	int ij,kl,ijkl,ik,jl,ikjl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			F(i,j) = scf.H_core(i,j);
//			cout << "before: " << F(i,j) << endl;
			for (int k = 0; k < 7; k++) {
				for (int l = 0; l < 7; l++) {
					ij = scf.index(i,j);
					kl = scf.index(k,l);
					ijkl = scf.index(ij,kl);
					ik = scf.index(i, k);
					jl = scf.index(j, l);
					ikjl = scf.index(ik, jl);
					F(i,j) += scf.D(k,l) * (2.0*scf.ERI(ijkl) - scf.ERI(ikjl));
				}
			}
			//cout << "after: " << F(i,j) << endl;
		}
	}
}

void	scf_procedure( HFSCF& scf, double delta1, double delta2 )
{
	Matrix	prev_D,F(7,7),D(7,7);
	double	prev_E,delta_E,rms;
	int	iter = 0;
	printf("Iter   \tE(elec)   \tE(tot)   \tDelta(E)   \tRMS(D)\n");
	printf("%d  %12.12f  %12.12f  \n",iter,scf.scf_E-scf.E_nuc,scf.scf_E);
	while (true) {
		iter++;
		prev_E = scf.scf_E;
		prev_D = scf.D;
		new_fock(scf, F);
		scf.DensityMatrix( F );
		scf.SCF_E( F );
		delta_E = scf.scf_E - prev_E;
		rms = sqrt(((scf.D-prev_D)*(scf.D-prev_D)).sum());
		printf("%d  %12.12f  %12.12f  %12.12f  %12.12f\n",iter,scf.scf_E-scf.E_nuc,scf.scf_E,delta_E,rms);
		if (delta_E<delta1 && rms<delta2) {
			break;
		}
	}
}

int main( int argc, char *argv[] ) {
	if (argc == 6) {
		HFSCF	scf( argv[1], argv[2], argv[3], argv[4], argv[5] );
		scf.compute_core_hamiltonian();
		scf.ortogon_S();
		scf.DensityMatrix( scf.H_core );
		scf.SCF_E( scf.H_core );
//		cout << "H_core: \n" << scf.H_core << endl;
		scf_procedure(scf, 1e-10, 1e-10 );
	}
}
