#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

int main(int, char**)
{
  cout.precision(3);
  cout << Matrix<double, 3, 4>::Identity() << endl;

  return 0;
}
