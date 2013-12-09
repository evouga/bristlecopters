#include <iostream>
#include <cmath>

using namespace std;

const double PI = 3.1415926535898;

int main()
{
	int numverts = 120;
	cout << numverts << " 2 0 0" << endl;
	for(int i=0; i<numverts; i++)
	{
		double x = cos( 2*PI*i/numverts );
		double y = sin( 2*PI*i/numverts );
		cout << i << " " << x << " " << y << endl;
	}
	cout << numverts << " 0" << endl;
	for(int i=0; i<numverts; i++)
		cout << i << " " << i << " " << (i+1)%numverts << endl;
	cout << "0" << endl;
}
