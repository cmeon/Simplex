#include <iostream>
#include <Eigen/Core>
#include "Simplex.h"
#include <Eigen/LU>
using namespace Eigen;
using namespace std;

typedef VectorXd ColVectorXd;

int main()
{
  int M=3; // number of constrains
  int N=3; // number of original variables
  MatrixXd A(M, N);

  double x[] = {3,1,3,6,2,2,3,3,0};
  double b[] = {22,14,14};
  double c[] = {1,4,5};
  // A.resize(M, N);
  A = Map<MatrixXd>(x, M, N);
  cout << A << endl;

  Simplex problem (M, N, x, c, b);
  problem.initTableau();
cout <<  problem.getTableau() << endl;
  problem.nextIteration();
    cout << problem.getTableau() << endl << endl;
   cout <<    problem.getCBt()<< endl;
   cout << problem.getB()<< endl;
   cout << problem.isOptimal()<< endl;
   cout << "c " << problem.getPivotCol()<< endl;
   cout << "r " << problem.getPivotRow()<< endl;
    
  problem.nextIteration();
    cout << problem.getTableau() << endl << endl;
   cout <<    problem.getCBt()<< endl;
   cout << "B= "<< problem.getB()<< endl;
   cout << problem.isOptimal()<< endl;
   cout << "c " << problem.getPivotCol()<< endl;
   cout << "r " << problem.getPivotRow()<< endl;
    
    
  problem.nextIteration();
    cout << problem.getTableau() << endl << endl;
   cout <<    problem.getCBt()<< endl;
   cout << problem.getB()<< endl;
   cout << problem.isOptimal()<< endl;
   cout << "c " << problem.getPivotCol()<< endl;
   cout << "r " << problem.getPivotRow()<< endl;
  return 0;
}
