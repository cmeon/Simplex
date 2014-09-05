// Simplex class

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

using namespace Eigen;

typedef VectorXd ColVectorXd; 
const   bool minimization = true;
const   bool maximization = false;

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

  bool natureOfProblem; // true if minimization, false if maximization

  int *xB;             //[M];
  RowVectorXd cBt;     //(M);
  MatrixXd B;          //(M, M);

  MatrixXd tableau;    //(M+2, M+N+2);

 public:
  Simplex(long int, long int, double[], double[], double[], int[], bool);
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
  long int getN();
  long int getM();
  long int getTableauRows();
  long int getTableauCols();
};

#endif // ifndef SIMPLEX_H
