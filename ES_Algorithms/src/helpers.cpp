
#include "HF_SCF.hpp"
vector <int>	ioff;

void	ioff_matrix( void )
{
	ioff.resize(BIGNUM);
	ioff[0] = 0;
	for (int i = 1; i < BIGNUM; i++) {
		ioff[i] = ioff[i-1] + i;
	}
}

int	index( int i, int j ) 
{
	int ij = i>j ? ioff[i]+j : ioff[j]+i;
	return ij;
}

void	parse_1e_integrals( char* file, Matrix& mtr )
{
	FILE	*f;
	int x,y;
	double val;
	f = fopen(file, "r");
	if (!f) {
		cout << "Unable to open: " << file << endl;
		exit(1);
	}
	while (fscanf(f, "%d %d %lf", &x,&y,&val) != EOF) {
		mtr(x-1,y-1) = val;
		mtr(y-1,x-1) = val;
	}
}

void	parse_eri( char *file, Array& ERI )
{
	FILE	*f;
	f = fopen(file, "r");
	int i,j,k,l,ij,kl,ijkl;
	double	val;
	if (!f) {
		cout << "Unable to open: " << file << endl;
		exit(1);
	}
	while (fscanf(f, "%d %d %d %d %lf", &i,&j,&k,&l,&val) != EOF) {
		ij = index(i-1, j-1);
		kl = index(k-1, l-1);
		ijkl = index(ij, kl);
		ERI(ijkl) = val;
	}
}

void	zero_out( Matrix& m )
{
	for (int i = 0; i < m.size()/7; i++) {
		for (int j = 0; j < m.size()/7; j++) {
			if (abs(m(i,j)) < 1e-10){
				m(i,j) = 0.00000000;
			}
		}

	}
}

void	mmult(Matrix& A, int a, Matrix& B, int b, Matrix& C, int l, int m, int n)
{
	int i,j,k;
	if (!a && !b) {
		for (i=0;i<l;i++) {
			for (j=0;j<m;j++) {
				for (k=0;k<n;k++) {
					C(i,j) += A(i,j) * B(k,j);
				}
			}
		}
	}
	else if (!a && b) {
		for (i=0;i<l;i++) {
			for (j=0;j<m;j++) {
				for (k=0;k<n;k++) {
					C(i,j) += A(i,k) * B(j,k);
				}
			}
		}
	}
	else if (a && !b) {
		for (i=0;i<l;i++) {
			for (j=0;j<m;j++) {
				for (k=0;k<n;k++) {
					C(i,j) += A(k,i) * B(k,j);
				}
			}
		}
	}
	if (a && b) {
		for (i=0;i<l;i++) {
			for (j=0;j<m;j++) {
				for (k=0;k<n;k++) {
					C(i,j) += A(k,i) * B(j,k);
				}
			}
		}
	}
}

void	zeros(Matrix& X)
{
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			X(i,j) = X(j,i) = 0.0;
		}

	}
}
