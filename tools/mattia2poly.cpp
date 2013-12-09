#include <iostream>
#include <vector>

using namespace std;

int main()
{
	vector<double> x, y;
	do
	{
		double curx, cury;
		cin >> curx >> cury;
		if(cin)
		{
			x.push_back(curx);
			y.push_back(cury);
		}
	}
	while(cin);

	int numpts = x.size();
	cout << numpts << " 2 0 0" << endl;
	for(int i=0; i<numpts; i++)
		cout << i << " " << x[i] << " " << y[i] << endl;
	cout << numpts << " 0" << endl;
	for(int i=0; i<numpts; i++)
		cout << i << " " << i << " " << (i+1)%numpts << endl;
	cout << "0" << endl;
	return 0;
}
