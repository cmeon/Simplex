// Simplex class

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

using namespace Eigen;

typedef VectorXd ColVectorXd; 

class Simplex
{
 private:
  long int M;          // number of constrains
  long int N;          // number of original variables

  RowVectorXd cT;      //(N);
  MatrixXd A;          //(M, N);
  ColVectorXd b;       //(M);

  long int pivotRow;
  long int pivotCol;

  int *xB;             //[M];
  RowVectorXd cBt;     //(M);
  MatrixXd B;          //(M, M);

  MatrixXd tableau;    //(M+2, M+N+2);

 public:
  Simplex(long int, long int, double[], double[], double[], int[]);
  ~Simplex();
  void setPivot();
  bool nextIteration();
  void initTableau();
  MatrixXd getTableau();
  int* getXB();
  MatrixXd getCBt();
  MatrixXd getB();
  bool isOptimal();
  int getPivotCol();
  int getPivotRow();
  bool hasMultipleOptimalSolutions();
  bool hasUnboundedSolutions();
};

#endif // ifndef SIMPLEX_H
