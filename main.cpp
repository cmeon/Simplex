#include <iostream>
#include "Simplex.h"
using namespace std;

typedef VectorXd ColVectorXd;

int main()
{
  int M=3; // number of constrains
  int N=3; // number of original variables

  double A[] = {3,1,3,6,2,2,3,3,0};
  double b[] = {22,14,14};
  double c[] = {1,4,5};
  int    x[] = {1,1,1};

  Simplex problem (M, N, A, c, b, x);

  problem.initTableau();
  cout << problem.getTableau() << endl << endl;
  while (problem.nextIteration()) {
    cout << problem.getTableau() << endl << endl;
    cout <<    problem.getCBt()<< endl;
    cout << "B= "<< problem.getB()<< endl;
    cout << problem.isOptimal()<< endl;
    cout << "c " << problem.getPivotCol()<< endl;
    cout << "r " << problem.getPivotRow()<< endl;

    if (problem.hasMultipleOptimalSolutions()) {
      break;
    }

    if (problem.hasUnboundedSolutions()) {
      break;
    }
  }
  return 0;
}
