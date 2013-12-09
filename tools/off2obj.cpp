#include <iostream>
#include <string>
using namespace std;

int main()
{
	string header;
	cin >> header;
	int numv, numf, dummy;
	cin >> numv >> numf >> dummy;
	for(int i=0; i<numv; i++)
	{
		double x,y,z;
		cin >> x >> y >> z;
		cout << "v " << x << " " << y << " " << z << endl;
	}
	for(int i=0; i<numf; i++)
	{
		int v1, v2, v3;
		cin >> dummy >> v1 >> v2 >> v3;		
		cout << "f " << 1+v1 << " " << 1+v2 << " " << 1+v3 << endl;
	}
}
