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

#include "pseudo_inverse.hpp"
#include "print.hpp"
#include <Eigen/LU>
#include <Eigen/SVD>
#include <iostream>


namespace kinematic_elastic {
  
  
  void pseudo_inverse_nonsingular (Matrix const & mx,
				   Matrix & inv)
  {
    Matrix mmt = mx * mx.transpose();
    inv = mx.transpose() * mmt.inverse();
  }
  
  
  /*
    Cobbled together from
    http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
    and
    http://eigen.tuxfamily.org/index.php?title=FAQ
  */
  void pseudo_inverse_moore_penrose (Matrix const & mx,
				     Matrix & inv)
  {
    Eigen::JacobiSVD<Matrix> svd(mx, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (0 == svd.nonzeroSingularValues()) {
      inv = Matrix::Zero(mx.cols(), mx.rows());
      return;
    }
    
    double thresh(numeric_limits<double>::epsilon() * svd.singularValues()[0]);
    if (mx.rows() > mx.cols()) {
      thresh *= mx.rows();
    }
    else {
      thresh *= mx.cols();
    }
    
    Vector rec(Vector::Zero(svd.singularValues().size()));
    for (ssize_t ii(0); ii < svd.nonzeroSingularValues(); ++ii) {
      if (svd.singularValues()[ii] <= thresh) {
	break;
      }
      rec[ii] = 1.0 / svd.singularValues()[ii];
    }
    
    inv = svd.matrixV() * rec.asDiagonal() * svd.matrixU().transpose();
  }
  
  
  void pseudo_inverse_damped (Matrix const & mx,
			      double lambda,
			      Matrix & inv)
  {
    if (lambda <= 0.0) {
      pseudo_inverse_moore_penrose (mx, inv);
      return;
    }
    
    Eigen::JacobiSVD<Matrix> svd(mx, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (0 == svd.nonzeroSingularValues()) {
      inv = Matrix::Zero(mx.cols(), mx.rows());
      return;
    }
    
    lambda *= lambda;
    
    inv
      = (svd.singularValues()[0] / (pow(svd.singularValues()[0], 2.0) + lambda))
      * svd.matrixV().col(0)
      * svd.matrixU().col(0).transpose();
    for (ssize_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
      inv
	+= (svd.singularValues()[ii] / (pow(svd.singularValues()[ii], 2.0) + lambda))
	* svd.matrixV().col(ii)
	* svd.matrixU().col(ii).transpose();
    }
  }
  
  
  void pseudo_inverse_baerlocher (Matrix const & mx,
				  double d_damp,
				  Matrix & inv,
				  Matrix & delta_projector,
				  Matrix * dbgU,
				  Vector * dbgsigma,
				  Matrix * dbgV,
				  Vector * dbgdamping)
  {
    Eigen::JacobiSVD<Matrix> svd(mx, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (0 == svd.nonzeroSingularValues()) {
      inv = Matrix::Zero(mx.cols(), mx.rows());
      delta_projector = Matrix::Zero(mx.cols(), mx.cols());
      return;
    }
    
    if (dbgU) {
      *dbgU = svd.matrixU();
    }
    if (dbgsigma) {
      *dbgsigma = svd.singularValues();
    }
    if (dbgV) {
      *dbgV = svd.matrixV();
    }
    
    typedef Eigen::JacobiSVD<Matrix>::Index index_t;
    
    double const sigma_min(svd.singularValues()[svd.nonzeroSingularValues() - 1]);
    if (sigma_min >= d_damp) {
      
      // no need for damping, use straight Moore-Penrose pseudo-inverse
      inv
	= (1.0 / svd.singularValues()[0])
	* svd.matrixV().col(0)
	* svd.matrixU().col(0).transpose();
      for (index_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
	inv
	  += (1.0 / svd.singularValues()[ii])
	  * svd.matrixV().col(ii)
	  * svd.matrixU().col(ii).transpose();
      }
      
      if (dbgdamping) {
	*dbgdamping = Vector::Zero(svd.nonzeroSingularValues() + 1);
	for (index_t ii(0); ii < svd.nonzeroSingularValues(); ++ii) {
	  (*dbgdamping)[ii] =
	    1.0 / svd.singularValues()[ii];
	}
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
      
      inv
	= (svd.singularValues()[0] / (pow(svd.singularValues()[0], 2.0) + lsquare))
	* svd.matrixV().col(0)
	* svd.matrixU().col(0).transpose();
      for (index_t ii(1); ii < svd.nonzeroSingularValues(); ++ii) {
	inv
	  += (svd.singularValues()[ii] / (pow(svd.singularValues()[ii], 2.0) + lsquare))
	  * svd.matrixV().col(ii)
	  * svd.matrixU().col(ii).transpose();
      }
      
      if (dbgdamping) {
	*dbgdamping = Vector::Zero(svd.nonzeroSingularValues() + 1);
	for (index_t ii(0); ii < svd.nonzeroSingularValues(); ++ii) {
	  (*dbgdamping)[ii] =
	    svd.singularValues()[ii] / (pow(svd.singularValues()[ii], 2.0) + lsquare);
	}
	(*dbgdamping)[svd.nonzeroSingularValues()] = sqrt(lsquare);
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
  
}
