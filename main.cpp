#include <iostream>
#include "Simplex.h"
#include <limits>
using namespace std;

typedef VectorXd ColVectorXd;

int main()
{
  int M=4; // number of constrains
  int N=2; // number of original variables

  double A[] = {1,2,-1,0,2,1,1,1};
  double b[] = {6,8,1,2};
  double c[] = {3,2};
  int    x[] = {1,1};
  bool min   = false;

  Simplex problem (M, N, A, c, b, x, min);

  problem.initTableau();
  cout << problem.getTableau() << endl << endl;
  do {
    problem.nextIteration();
    cout << problem.getTableau() << endl << endl;

    if (problem.hasMultipleOptimalSolutions()) {
      cout << "multiple" << endl;
      break;
    }

    if (problem.hasUnboundedSolutions() && !problem.isOptimal()) {
      cout << "unbounded" << endl;
      break;
    }
  }   while (!problem.isOptimal()) ;
  //using a built-in function to display the machine-epsilon given the data type
  std::cout << "The machine precision for double is : " << std::numeric_limits<double>::epsilon() << std::endl;
  
  return 0;
}
