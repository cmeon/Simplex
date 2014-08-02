// Simplex class

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

using namespace Eigen;

typedef VectorXd ColVectorXd; 

class Simplex
{
 private:
  long int M; // number of constrains
  long int N; // number of original variables

  RowVectorXd cT; //(N);
  MatrixXd A; //(M, N);
  ColVectorXd b; //(M);
  struct {
    int row, col;
  } pivot;
  long int pivotRow;
  long int pivotCol;
  int *xB; // [M];
  RowVectorXd cBt; //(M);
  MatrixXd B; //(M, M);

  MatrixXd tableau; //(M+2, M+N+2);

 public:
  // constructor
  Simplex(long int, long int, double[], double[], double[]);
  ~Simplex();
  void setPivot();

 public:
  bool nextIteration();
  void initTableau();
  MatrixXd getTableau();
  int* getXB();
  MatrixXd getCBt();
  MatrixXd getB();
  bool isOptimal();
  int getPivotCol();
  int getPivotRow();
};
