#include "MP2.hpp"
#include "HF_SCF.hpp"
#include <chrono>

MP2::MP2( Array& eri, Matrix& c, Array& E_mos ) : AO_ERI(eri), C(c), E_MO(E_mos)
{
	MO_ERI.resize((7*(7+1)/2)*((7*(7+1)/2)+1)/2);
	E_mp2 = 0.0;
}

MP2::~MP2( void )
{}

void	MP2::Noddy_MO_Transformation( void )
{
	int i,j,k,l,ijkl;
	int p,q,r,s,pq,rs,pqrs;
	for (i=0,ijkl=0;i<7;i++) {
		for (j=0;j<=i;j++) {
			for (k=0;k<=i;k++) {
				for (l=0;l<(i==k?j:k);l++,ijkl++) {
					for (int p = 0; p < 7; p++) {
						for (int q = 0; q < 7; q++) {
							pq = index(p,q);
							for (int r = 0; r < 7; r++) {
								for (int s = 0; s < 7; s++) {
									rs = index(r,s);
									pqrs = index(pq,rs);
									MO_ERI(ijkl) += C(p,i)*C(q,j)*C(r,k)*C(s,l)*AO_ERI(pqrs);
								}
							}
						}
					}

				}
			}
		}
	}
}

void	MP2::Smarter_MO_Transformation( void )
{
	Matrix	X(7,7),Y(7,7),TMP(7*(7+1)/2,7*(7+1)/2);
	int i,j,k,l,ij,kl,ijkl,klij;
	for (i=0,ij=0;i<7;i++) {
		for (j=0;j<=i;j++,ij++) {
			for (k=0,kl=0;k<7;k++) {
				for (l=0;l<=k;l++,kl++) {
					ijkl = index(ij, kl);
					X(k,l) = X(l,k) = AO_ERI(ijkl);
				}
			}
			zeros(Y);
			mmult(C, 1, X, 0, Y, 7, 7, 7);
			zeros(X);
			mmult(Y, 0, C, 0, X, 7, 7, 7);
			for (k=0,kl=0;k<7;k++) {
				for (l=0;l<=k;l++,kl++) {
					TMP(kl,ij) = X(k,l);
					}
				}
		}
	}
	for (k=0,kl=0;k<7;k++) {
		for (l=0;l<=k;l++,kl++) {
			zeros(X);
			zeros(Y);
			for (i=0,ij=0;i<7;i++) {
				for (j=0;j<=i;j++,ij++) {
					X(i,j) = X(j,i) = TMP(kl,ij);
				}
			}
			zeros(Y);
			mmult(C, 1, X, 0, Y, 7, 7, 7);
			zeros(X);
			mmult(Y, 0, C, 0, X, 7, 7, 7);
			for (i=0,ij=0;i<7;i++) {
				for (j=0;j<=i;j++,ij++) {
					klij = index(kl, ij);
					MO_ERI(klij) = X(i,j);
					}
				}
		}
	}
}

void	MP2::computeE( void )
{
	auto start = std::chrono::high_resolution_clock::now();
	Smarter_MO_Transformation();
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
	cout << "Noddy took: " << duration.count() << " ms" << endl;
	int i,j,a,b,ia,jb,ib,ja,iajb,ibja;
	for (i=0;i<5;i++) {
		for (a=5;a<7;a++) {
			ia = index(i, a);
			for (j=0;j<5;j++) {
				ja = index(j, a);
				for (b=5;b<7;b++) {
					jb = index(j, b);
					ib = index(i, b);
					iajb = index(ia, jb);
					ibja = index(ib, ja);
					E_mp2 += MO_ERI(iajb)*(2*MO_ERI(iajb)-MO_ERI(ibja))/(E_MO(i)+E_MO(j)-E_MO(a)-E_MO(b));
				}
			}
		}
	}
}

double&	MP2::E_MP2( void )
{
	return E_mp2;
}
