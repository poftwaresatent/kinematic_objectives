// g++ -I/usr/include/eigen3 -Wall -o tstswap tstswap.cpp

#include <iostream>
#include <Eigen/Core>
#include <err.h>

using namespace std;
using namespace Eigen;

int main(int argc, char ** argv)
{
  MatrixXd aa, bb, cc, dd;
  
  aa = MatrixXd::Random(17, 42);
  bb = MatrixXd::Random(33, 22);
  cc = aa;
  dd = bb;
  
  aa.swap(bb);
  
  if (aa.rows() != dd.rows()) {
    errx(EXIT_FAILURE, "aa.rows() == %zd != dd.rows() == %zd", aa.rows(), dd.rows());
  }
  if (aa.cols() != dd.cols()) {
    errx(EXIT_FAILURE, "aa.cols() == %zd != dd.cols() == %zd", aa.cols(), dd.cols());
  }
  for (ssize_t ii(0); ii < aa.rows(); ++ii) {
    for (ssize_t jj(0); jj < aa.cols(); ++jj) {
      if (aa(ii, jj) != dd(ii, jj)) {
	errx(EXIT_FAILURE, "aa(%zd, %zd) != dd(%zd, %zd)", ii, jj, ii, jj);
      }
    }
  }
  
  if (bb.rows() != cc.rows()) {
    errx(EXIT_FAILURE, "bb.rows() == %zd != cc.rows() == %zd", bb.rows(), cc.rows());
  }
  if (bb.cols() != cc.cols()) {
    errx(EXIT_FAILURE, "bb.cols() == %zd != cc.cols() == %zd", bb.cols(), cc.cols());
  }
  for (ssize_t ii(0); ii < bb.rows(); ++ii) {
    for (ssize_t jj(0); jj < bb.cols(); ++jj) {
      if (bb(ii, jj) != cc(ii, jj)) {
	errx(EXIT_FAILURE, "bb(%zd, %zd) != cc(%zd, %zd)", ii, jj, ii, jj);
      }
    }
  }
  
  cout << "ok then, it seems to safe to use swap\n";
}
