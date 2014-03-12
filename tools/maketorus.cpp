#include <iostream>
#include <cmath>

using namespace std;

const double PI = 3.1415926535898;

int main()
{
	int numverts = 500;
	double holeradius = 0.1;
	int numintverts = int(holeradius*numverts);
	cout << numverts+numintverts << " 2 0 0" << endl;
	for(int i=0; i<numverts; i++)
	{
		double x = cos( 2*PI*i/numverts );
		double y = sin( 2*PI*i/numverts );
		cout << i << " " << x << " " << y << endl;
	}
	for(int i=0; i<numintverts; i++)
	{
		double x = holeradius*cos( 2*PI*i/numintverts );
		double y = holeradius*sin( 2*PI*i/numintverts );
		cout << i+numverts << " " << x << " " << y << endl;
	}
	cout << numverts+numintverts << " 0" << endl;
	for(int i=0; i<numverts; i++)
		cout << i << " " << i << " " << (i+1)%numverts << endl;
	for(int i=0; i<numintverts; i++)
		cout << i+numverts << " " << i+numverts << " " << numverts+((i+1)%numintverts) << endl;
	cout << "1" << endl;
	cout << "0 0 0" << endl;
}
