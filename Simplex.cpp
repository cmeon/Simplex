#include "Simplex.h"
#include <iostream>
#include <cmath>

/*
  Simplex problem constructor

  This takes the one dimensional arrays:
    A - the LHS variable coefficient matrix of the constraints,
    c - the row vector for the coeffients of the variables,
    b - the column vector for the RHS

  It generates a simplex tableau by mapping the values accordingly
  and calling initTableau to create the simplex tableau.

  Array x stores the nature of the variables such that a value of:
    -1 - corresponds to x>=0
     0 - means x is unrestricted
    +1 - means x<=0

  m is number of constraints
  n is number of variables

*/

Simplex::Simplex(long int m, long int n, double A[], double c[], double b[], int x[], bool natureOfProblem)
{

  this->natureOfProblem = natureOfProblem;

  if (natureOfProblem == maximization) {
    int colSize = n; // size of the matrix A column
    int rowSize = m; // size of the matrix A row
    
    // to create a modified version of A and c, we calculate the size
    // of the new modified A and c. Only the column size will change
    // depending on the nature of the variables
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
    
    // laundry work!
    delete[] modifiedA;
    delete[] modifiedC;
  }

  if (natureOfProblem == minimization) {
    // create a temporary matrix
    MatrixXd temp(m+1, n+1);

    double* modifiedA = new double[m*n];
    double* modifiedC = new double[n];
    
    int col = 0;
    
    int colSize = n; // size of the matrix A column
    int rowSize = m; // size of the matrix A row
    
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
    
    temp.block(0, 0, m, n) = Map<MatrixXd>(modifiedA, m, n);
    temp.block(m, 0, 1, n) = Map<RowVectorXd>(modifiedC, n);
    temp.block(0, n, m, 1) = Map<ColVectorXd>(b, m);
    temp(m, n)               = 0;
    MatrixXd tempTransposed = temp.transpose();
    
    double* tempA = new double[m*n];
    for (int i=0; i<n; i++) {
      for (int j=0; j<m; j++) {
	tempA[j*n+i] = tempTransposed(i,j);
      }
      std::cout << '\n';
    }
    
    double* tempb = new double[n];
    for (int i=0; i<n; i++) {
      tempb[i] = tempTransposed(i, m);
    }
    
    double* tempc = new double[m];
    for (int i=0; i<m; i++) {
      tempc[i] = tempTransposed(n, i);
    }
    
    this->M = n;
    this->N = m;
    
    this->A  = Map<MatrixXd>(tempA, n, m);
    this->cT = Map<RowVectorXd>(tempc, m);
    this->b  = Map<ColVectorXd>(tempb, n);
    
    // laundry work!
    delete[] modifiedA;
    delete[] modifiedC;
    
    delete[] tempA;
    delete[] tempb;
    delete[] tempc;
  }

  initTableau();
}



/*
  Generate initial tableau for basic feasible solution
*/

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

  // initialize xB which holds the basic variable column number
  xB = new int[M];
  cBt.resize(M);

  for (int i = 0; i<M; i++) {
    xB[i]  = i+N+1;
    cBt[i] = tableau.row(0)[xB[i]];
  }

  // creating the B matrix which is the initial basic matrix
  B = tableau.block(1, N+1, M, M);

  // set pivot
  setPivot();
}



/*
  Set the pivot location
*/

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



/*
  Calculate the next iteration

  It returns true if the there is a next iteration and false if the Simplex
  problem is solved
*/

bool Simplex::nextIteration()
{
  
  setPivot();

  xB[pivotRow-1]    = pivotCol;
  cBt[pivotRow-1]   = cT[pivotCol-1];
  B.col(pivotRow-1) = A.col(pivotCol-1);

  MatrixXd Binv = B.inverse();

  tableau.block(0, 1, 1, N)     = ( (cBt * Binv * A) - cT );
  tableau.block(1, 1, M, N)     = ( Binv * A );
  tableau.block(0, N+1, 1, M)   = ( cBt * Binv );
  tableau.block(1, N+1, M, M)   = ( Binv );
  tableau(0, M+N+1)             = ( cBt * Binv * b );
  tableau.block(1, M+N+1, M, 1) = ( Binv * b );

  if (isOptimal())
    return false;

  return true;
}



/*
  Get the Simplex tableau

  Returns a matrix object(from eigen implementaion) with the simplex tableau
*/

MatrixXd Simplex::getTableau()
{
  return tableau;
}



/*
  Get the values of xB

  Returns a copy of the basic variables
*/

int* Simplex::getXB()
{
  // copy the contents of 'xB' in new 'x' to protect memory from being exposed
  // in the outside environment
  int* x = new int[M];
  for (int i=0; i<M; i++)
    x[i] = xB[i];
  return x;
}



/*
  Get the current cBt matrix

  Returns a Matrix instance that contains the values of cBt
*/

MatrixXd Simplex::getCBt()
{
  return cBt;
}



/*
  Get the current B matrix

  Returns a Matrix instance with values B
*/

MatrixXd Simplex::getB()
{
  return B;
}



/*
  Optimality check

  Returns true if the matrix is optimal
*/

bool Simplex::isOptimal()
{
  bool result = true;
  RowVectorXd temp = tableau.block(0, 1, 1, N);
  int *nonBasic = new int[N];
  RowVectorXd zRow = this->tableau.block(0, 1, 1, N);
  for (int i=0; i<M+N; i++) {
    nonBasic[i] = 0;
  }

  for (int i=0; i<M; i++) {
    nonBasic[xB[i]-1] = 1;
  }

  for (int i=0; i<N; i++) {
    if (nonBasic[i] == 0) {
      if (zRow(i) > 0) {
	result &= true;
      } else {
	result &= false;
      }
    }
    std::cout << nonBasic[i] << ", " << std::endl;
    if (nonBasic[i] == 1) {
      if (zRow(i) == 0) {
	result &= true;
      } else {
	result &= false;
      }
    }
  }

  /*
  for (int i =0; i<N; i++) {
    if (!(std::signbit(temp[i]) == 0)) {
      return false;
    }
  }
  
  return true;
  */
  return result;
}



/*
  Get the pivot column

  Returns the value of pivot column corresponding to the simplex tableaux
*/

int Simplex::getPivotCol()
{
  return pivotCol;
}



/*
  Get the pivot row

  Returns the value of pivot row corresponding to the simplex tableaux
*/

int Simplex::getPivotRow()
{
  return pivotRow;
}



/*
  Check for case of unbounded solutions

  The situation occurs when the coefficient for the denominator of 
  the intercept ratios are either zero or negative.
  Returns true or false true indicating the unbounded solution
*/

bool Simplex::hasUnboundedSolutions()
{
  ColVectorXd pivotColumn = tableau.block(1, pivotCol-1, M, 1);

  bool result = true;

  if ((pivotColumn.array() > (-std::numeric_limits<double>::epsilon())).all()) {
    result = false;
  }

  return result;
}



/*
  Check for case of Multiple optimal solution / Alternative optimal solutions

  When zero appears in the column of a nonbasic variable in the
  z-row of the optimal tableau
  iReturns true or false true indicating the unbounded solution
*/

bool Simplex::hasMultipleOptimalSolutions() {
  int *nonBasic = new int[M+N];

  RowVectorXd zRow = this->tableau.block(0, 1, 1, N+M);
  for (int i=0; i<M+N; i++) {
    nonBasic[i] = 0;
  }

  for (int i=0; i<M; i++) {
    nonBasic[xB[i]-1] = 1;
  }

  for (int i=0; i<M+N; i++) {
    if (nonBasic[i] == 0) {
      if (zRow(i) == 0) {
	      return true;
      }
    }  
  }
  
  delete [] nonBasic;
  return false;
  
}



/*
  gets N
*/
long int Simplex::getN() {
  return this->N;
}



/*
  gets M
*/
long int Simplex::getM() {
  return this->M;
}



/*
  get tableau size
*/
long int Simplex::getTableuCols() {
  return tableau.cols();
}



/*
  get tableau size
*/
long int Simplex::getTableuRows() {
  return tableau.rows();
}



/*
  Destructor
*/

Simplex::~Simplex()
{
  delete[] xB;
}

