
#include "HF_SCF.hpp"
#include <cmath>
#include <cstdio>

void	new_fock( HFSCF& scf )
{
	int ij,kl,ijkl,ik,jl,ikjl;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			scf.F(i,j) = scf.H_core(i,j);
			for (int k = 0; k < 7; k++) {
				for (int l = 0; l < 7; l++) {
					ij = scf.index(i,j);
					kl = scf.index(k,l);
					ijkl = scf.index(ij,kl);
					ik = scf.index(i, k);
					jl = scf.index(j, l);
					ikjl = scf.index(ik, jl);
					scf.F(i,j) += scf.D(k,l) * (2.0*scf.ERI(ijkl) - scf.ERI(ikjl));
				}
			}
		}
	}
}

void	scf_procedure( HFSCF& scf, double delta1, double delta2 )
{
	Matrix	prev_D;
	double	prev_E,delta_E,rms;
	int	iter = 0;
	printf("Iter   \tE(elec)   \tE(tot)   \tDelta(E)   \tRMS(D)\n");
	while (true) {
		prev_E = scf.scf_E;
		prev_D = scf.D;
		new_fock(scf);
		scf.DensityMatrix( scf.F );
		scf.SCF_E( scf.F );
		delta_E = scf.scf_E - prev_E;
		rms = sqrt(((scf.D-prev_D)*(scf.D-prev_D)).sum());
		printf("%d  %12.12f  %12.12f  %12.12f  %12.12f\n",iter,scf.scf_E-scf.E_nuc,scf.scf_E,delta_E,rms);
		if (delta_E<delta1 && rms<delta2) {
			break;
		}
		iter++;
	}
}

int main( int argc, char *argv[] ) {
	if (argc == 6) {
		HFSCF	scf( argv[1], argv[2], argv[3], argv[4], argv[5] );
		scf.compute_core_hamiltonian();
		scf.ortogon_S();
		scf.DensityMatrix( scf.H_core );
		scf.SCF_E( scf.H_core );
		cout << "HAHA\n";
		scf_procedure(scf, 1e-10, 1e-10 );
	}
}
