
#include "HF_SCF.hpp"
#include "MP2.hpp"
#include <cstdio>

void	MO_Fock( Matrix& AOF, Matrix& C )
{
	Matrix	MF(7,7);
	for (int p = 0; p < 7; p++) {
		for (int q = 0; q < 7; q++) {
			MF(p,q) = 0.0;
			for (int r = 0; r < 7; r++) {
				for (int s = 0; s < 7; s++) {
					MF(p,q) += C(r,q) * C(s,p) * AOF(r,s);
				}
			}
	}
	}
	cout << endl << MF << endl;
}

int main( int argc, char *argv[] ) {
	if (argc == 6) {
		ioff_matrix();
		SCF	scf( argv[1], argv[2], argv[3], argv[4], argv[5] );
		scf.initialize();
		scf.runSCF(1e-10, 1e-10);
		MO_Fock(scf.AOFock(), scf.CoefficientMatrix());
		MP2	mp2(scf.AOERI(), scf.CoefficientMatrix(), scf.MO_Energies());
		mp2.computeE();
		printf("MP2 energy: %12.12f\n", mp2.E_MP2());
		printf("Total energy: %12.12f\n", scf.E() + mp2.E_MP2());
	}
}
