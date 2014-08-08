#include "Simplex.h"
#include <iostream>
#include <cmath>

// Simplex problem constructor
Simplex::Simplex(long int m, long int n, double A[], double c[], double b[], int x[])
{
  int colSize = n;
  int rowSize = m;
  for (int i = 0; i<m; i++) {
    if (x[i] == 0) n++;
  }
  this->M = m;
  this->N = n;
  double* modifiedA = new double[m*n];
  double* modifiedC = new double[n];

  int col = 0;

  for (int i=0; i<colSize; i++) {
    for (int j=0; j<rowSize; j++) {
      int idx = i*m+j;
      int idxMA = col*m+j;

      if (x[i] ==-1) {
	modifiedA[idxMA]   = -(A[idx]);
	modifiedC[col] = -(c[i]);
      }
      if (x[i] == 0) {
	modifiedA[idxMA]   = +(A[idx]);
	modifiedA[idxMA+m] = -(A[idx]);
	modifiedC[col] = +(c[i]);
	modifiedC[col+1] = -(c[i]);
      }
      if (x[i] == 1) {
	modifiedA[idxMA]   = +(A[idx]);
	modifiedC[col] = +(c[i]);
      }
    }

    if (x[i] == 0) col+=2;
    else           col++;

  }

  this->A  = Map<MatrixXd>(modifiedA, m, n);
  this->cT = Map<RowVectorXd>(modifiedC, n);
  this->b  = Map<ColVectorXd>(b, m);

  delete[] modifiedA;
  delete[] modifiedC;

  initTableau();
}

// Generate initial tableau for basic feasible solution
void Simplex::initTableau()
{
  this->tableau.resize(M+1, M+N+2);
  this->tableau(0, 0)               = 1;
  this->tableau.block(0, 1, 1, N)   = (this->cT.array() * -1);
  this->tableau.block(0, N+1, 1, M) = RowVectorXd::Zero(M);
  this->tableau(0, N+1)             = 0;

  this->tableau.block(1, 0, M, 1)   = ColVectorXd::Zero(M);
  this->tableau.block(1, 1, M, N)   = this->A;
  this->tableau.block(1, N+1, M, M) = MatrixXd::Identity(M, M);

  this->tableau.block(1,M+N+1,M, 1) = this->b;

  // initialize xB
  xB = new int[M];
  cBt.resize(M);

  for (int i = 0; i<M; i++) {
    xB[i]  = i+N+1;
    cBt[i] = tableau.row(0)[xB[i]];
  }

  B = tableau.block(1, N+1, M, M);
}

// Set the pivot location
void Simplex::setPivot()
{
  double min = tableau(0, 1);
  for (int i=1; i<=N; i++) {
    if (tableau(0, i) <= min) {
      min = tableau(0, i);
      pivotCol = i;
    }
  }
  ColVectorXd solution    = tableau.block(1, M+N+1,    M, 1);
  ColVectorXd pivotColumn = tableau.block(1, pivotCol, M, 1);

  ArrayXd ratio = solution.array() / pivotColumn.array();
  min = ratio(0);
  pivotRow = 0;
  for (int i=0; i<M; i++) {
    if (ratio(i) <= min) {
      min = ratio(i);
      pivotRow = i+1;
    }
  }
}

// update simplex xB list
bool Simplex::nextIteration()
{
  if (isOptimal())
    return false;

  setPivot();

  xB[pivotRow-1]    = pivotCol-1;
  cBt[pivotRow-1]   = cT[pivotCol-1];
  B.col(pivotRow-1) = A.col(pivotCol-1);

  MatrixXd Binv = B.inverse();

  tableau.block(0, 1, 1, N)     = ( (cBt * Binv * A) - cT );
  tableau.block(1, 1, M, N)     = ( Binv * A );
  tableau.block(0, N+1, 1, M)   = ( cBt * Binv );
  tableau.block(1, N+1, M, M)   = ( Binv );
  tableau(0, M+N+1)               = ( cBt * Binv * b );
  tableau.block(1, M+N+1, M, 1) = ( Binv * b );

  return true;
}


// get the Simplex tableau
MatrixXd Simplex::getTableau()
{
  return tableau;
}

// get the values of xB
int* Simplex::getXB()
{
  // copy the contents of 'xB' in new 'x' to protect memory from being exposed
  // in the outside environment
  int* x = new int[M];
  for (int i=0; i<M; i++)
    x[i] = xB[i];
  return x;
}

// Get the current cBt matrix
MatrixXd Simplex::getCBt()
{
  return cBt;
}

// Get the current B matrix
MatrixXd Simplex::getB()
{
  return B;
}

// Optimality check
bool Simplex::isOptimal()
{
  RowVectorXd temp = tableau.block(0, 1, 1, N);
  for (int i =0; i<N; i++) {
    if (!(std::signbit(temp[i]) == 0)) {
      return false;
    }
  }
  return true;
}

// Get the pivot column
int Simplex::getPivotCol()
{
  return pivotCol;
}

// Get the pivot row
int Simplex::getPivotRow()
{
  return pivotRow;
}

// Empty Destructor
Simplex::~Simplex()
{
  std::cout << xB[0] << std::endl;
  delete[] xB;
  std::cout << xB[0] << std::endl;
}
