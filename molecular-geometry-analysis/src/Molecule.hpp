

#include <map>
#include <string>
#include <vector>
using namespace std;

class	Molecule
{
	public:
		Molecule(string filename, int q);
		~Molecule();

	// NOTE: Attributes of a molecule
		int	natoms;
		int	charge;
		vector<int>	zvals;
		vector<vector<double> >	geom;
		string	point_group;

	// NOTE: Molecule members
		void	print_geom();
		void	rotate( double );
		void	translate( double , double, double );
		double	unit( int, int, int );
		double	bond( int, int );	
		double	angle( int, int, int );
		double	torsion( int, int, int, int );
		double	oop_angle( int, int, int, int );
};
