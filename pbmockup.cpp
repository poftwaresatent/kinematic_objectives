/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2013, Willow Garage, Inc.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Willow Garage nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Roland Philippsen */

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/StdVector>
//#include <vector>

#include <iostream>

#include <err.h>


namespace pbmockup {
  
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  
  using namespace std;
  
  
  struct task_s {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    task_s(size_t ndof, size_t ndim_, double b_max_)
      : b_max(b_max_),
	ndim(ndim_)
    {
      current = Vector::Zero(ndim_);
      desired = Vector::Zero(ndim_);
      Jacobian = Matrix::Zero(ndim_, ndof);
    }
    
    Vector current;
    Vector desired;
    Matrix Jacobian;
    double b_max;
    size_t ndim;
  };
  
  typedef vector<task_s> tasklist_t;
  
  
  struct system_s {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    system_s (size_t ndof_)
      : ndof(ndof_)
    {
      state = Vector::Zero(ndof_);
    }
    
    Vector state;
    size_t ndof;
  };
  
  
  void compute_svd_stuff(Matrix const & Jt,
			 double d_damp,
			 Matrix & Jt_inv_damped,
			 Matrix & delta_projector)
  {
    Eigen::JacobiSVD<Matrix> svd(Jt, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (0 == svd.nonzeroSingularValues()) {
      Jt_inv_damped = Matrix::Zero(Jt.cols(), Jt.rows());
      delta_projector = Matrix::Zero(Jt.cols(), Jt.cols());
      return;
    }
    
    typedef Eigen::JacobiSVD<Matrix>::Index index_t;
    
    double const sigma_min(svd.singularValues()[svd.nonzeroSingularValues() - 1]);
    if (sigma_min >= d_damp) {
      
      // no need for damping, use straight Moore-Penrose pseudo-inverse
      Jt_inv_damped
	= (1.0 / svd.singularValues()[0])
	* svd.matrixV().col(0)
	* svd.matrixU().col(0).transpose();
      for (index_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
	Jt_inv_damped
	  += (1.0 / svd.singularValues()[ii])
	  * svd.matrixV().col(ii)
	  * svd.matrixU().col(ii).transpose();
      }
    }
    else {
      
      double lsquare;
      if (sigma_min <= d_damp / 2.0) {
	lsquare = pow(d_damp / 2.0, 2.0);
      }
      else {
	lsquare = sigma_min * (d_damp - sigma_min);
      }
      
      Jt_inv_damped
	= (svd.singularValues()[0] / (pow(svd.singularValues()[0], 2.0) + lsquare))
	* svd.matrixV().col(0)
	* svd.matrixU().col(0).transpose();
      for (index_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
	Jt_inv_damped
	  += (svd.singularValues()[ii] / (pow(svd.singularValues()[ii], 2.0) + lsquare))
	  * svd.matrixV().col(ii)
	  * svd.matrixU().col(ii).transpose();
      }
    }
    
    // the following could be sped up because it produces a symmetric matrix
    delta_projector
      = svd.matrixV().col(0)
      * svd.matrixV().col(0).transpose();
    for (index_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
      delta_projector
	+= svd.matrixV().col(ii)
	* svd.matrixV().col(ii).transpose();
    }
  }
  
  
  Vector recursive_task_priority_algorithm (system_s const & system,
					    tasklist_t const & tasklist)
  {
    // initialization (without constraints for now...)
    Vector delta_q(Vector::Zero(system.ndof));
    Matrix projector(Matrix::Identity(system.ndof, system.ndof));
    
    // recursively add tasks, projected into the nullspace of
    // higher-priority tasks [Baerlocher 2001 Fig 5.12]
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      
      Vector dx = tasklist[ii].desired - tasklist[ii].current;
      double const d_damp(dx.norm() / tasklist[ii].b_max);
      ////dx /= d_damp;		// this is never explicitly stated in [Baerlocher, 2001]
      Vector dx_comp = dx - tasklist[ii].Jacobian * delta_q;
      
      Matrix Jtilda = tasklist[ii].Jacobian * projector;
      
      Matrix Jtilda_inv, delta_projector;
      compute_svd_stuff (Jtilda, d_damp, Jtilda_inv, delta_projector);
      
      delta_q = delta_q + Jtilda_inv * dx_comp;
      
      projector = projector - delta_projector;
    }
    
    // // add joint-space "criterion minimization term" (i.e. posture)
    
    // double const xi(-1.0e-4);	// must be negative
    // Vector const & posture_gradient(system.state); // Vector::Ones(ndof)); // some bias...
    
    // delta_q = delta_q + projector * xi * posture_gradient;
    
    return delta_q;
  }
  
}


using namespace pbmockup;


void dump (system_s const & system,
	   tasklist_t const & tasklist,
	   Vector const & dq)
{
  for (size_t ii(0); ii < system.ndof; ++ii) {
    cout << system.state[ii] << "  ";
  }
  for (size_t ii(0); ii < tasklist.size(); ++ii) {
    cout << "  ";
    for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
      cout << tasklist[ii].current[jj] << "  ";
    }
  }
  cout << "  ";
  for (size_t ii(0); ii < system.ndof; ++ii) {
    cout << dq[ii] << "  ";
  }
  cout << "\n";
}


void dbg (system_s const & system,
	  tasklist_t const & tasklist,
	  Vector const & dq)
{
  cout << "==================================================\n"
       << "state:";
  for (size_t ii(0); ii < system.ndof; ++ii) {
    cout << "\t" << system.state[ii];
  }
  for (size_t ii(0); ii < tasklist.size(); ++ii) {
    cout << "\ntask " << ii << "\n";
    cout << "  current:";
    for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
      cout << "\t" << tasklist[ii].current[jj];
    }
    cout << "\n  desired:";
    for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
      cout << "\t" << tasklist[ii].desired[jj];
    }
    cout << "\n  Jacobian:";	// hardcoded for 1xN matrices
    for (size_t jj(0); jj < system.ndof; ++jj) {
      cout << "\t" << tasklist[ii].Jacobian(0, jj);
    }
  }
  cout << "\ndelta_q:";
  for (size_t ii(0); ii < system.ndof; ++ii) {
    cout << "\t" << dq[ii];
  }
  cout << "\n";
}


int main (int argc, char ** argv)
{
  size_t const ndof(2);
  
  system_s system(ndof);
  system.state << 0.3, -0.2;
  
  tasklist_t tasklist;
  tasklist.push_back(task_s(ndof, 1, 1.0e-3));
  tasklist.push_back(task_s(ndof, 1, 1.0e-2 * M_PI / 180.0));
  
  tasklist[0].desired << 1.2;
  tasklist[1].desired << 35.0 * M_PI / 180.0;
  
  tasklist[1].Jacobian << 0.0, 1.0; // constant in this case
  
  for (;;) { ////size_t ii(0); ii < 10000; ++ii) {
    double const q0(system.state.coeff(0));
    double const q1(system.state.coeff(1));
    
    tasklist[0].current <<
      cos(q0) + cos(q0 + q1);
    tasklist[0].Jacobian <<
      -sin(q0) - sin(q0 + q1),
      -sin(q0 + q1);
    
    tasklist[1].current <<
      q1;
    
    Vector dq = recursive_task_priority_algorithm (system, tasklist);
    
    dump(system, tasklist, dq);
    
    system.state += dq;
    
  }
  
  return 0;
}
