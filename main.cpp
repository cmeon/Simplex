#include <iostream>
#include "Simplex.h"
using namespace std;

typedef VectorXd ColVectorXd;


int main()
{

  int M1 = 3; // number of constrains
  int N1 = 2; // number of original variables

  double A1[] = {4,4,4,3,1,-1};
  double b1[] = {12,8,8};
  double c1[] = {3,2};
  int    x1[] = {1,1};

  int M2 = 3;
  int N2 = 3;
  double A2[] = {3,1,3,6,2,2,3,3,0};
  double b2[] = {22,14,14};
  double c2[] = {1,4,5};
  int    x2[] = {1,1,1};
  
  //Simplex problem (M1, N1, A1, c1, b1, x1, false);
  Simplex problem (M2, N2, A2, c2, b2, x2, false);
  
  problem.initTableau();

  cout << problem.getTableau() << endl << endl;
  cout << "pivot =>" << problem.getPivotCol() << ", " << problem.getPivotRow() << endl;
  cout << "cbt   =>" << problem.getCBt()      << endl;
  cout << "B     =>" << problem.getB()        << endl;
  
  do {
    problem.nextIteration();
    
    cout << problem.getTableau() << endl << endl;
    cout << "pivot =>" << problem.getPivotCol() << ", " << problem.getPivotRow() << endl;
    cout << "cbt   =>" << problem.getCBt()      << endl;
    cout << "B     =>" << problem.getB()        << endl;
    
    if (problem.hasMultipleOptimalSolutions()) {
      cout << "multiple optimal solution.\n";
      //break;
    }
    
    if (problem.hasUnboundedSolutions()) {
      cout << "unbounded solution\n";
      break;
    }

    if (problem.isOptimal() && !problem.hasMultipleOptimalSolutions()) {
      cout << "optimal solution\n";
    }

  } while(!problem.isOptimal());
  return 0;
}

