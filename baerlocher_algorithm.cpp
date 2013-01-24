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

#include "baerlocher_algorithm.hpp"
#include "pseudo_inverse.hpp"
#include "print.hpp"
#include <Eigen/SVD>
#include <iostream>


namespace kinematic_elastic {
  
  
  Vector baerlocher_algorithm (Model const & model,
			       Vector const & state,
			       tasklist_t const & tasklist,
			       ostream * dbgos,
			       char const * dbgpre)
  {
    Matrix * dbgU(0);
    Matrix * dbgV(0);
    Vector * dbgsigma(0);
    Vector * dbgdamping(0);
    if (dbgos) {
      dbgU = new Matrix();
      dbgV = new Matrix();
      dbgsigma = new Vector();
      dbgdamping = new Vector();
    }
    
    // initialization (without constraints for now...)
    size_t const ndof(state.size());
    Vector delta_q(Vector::Zero(ndof));
    Matrix projector(Matrix::Identity(ndof, ndof));
    
    // recursively add tasks, projected into the nullspace of
    // higher-priority tasks [Baerlocher 2001 Fig 5.12]
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      
      Vector dx = tasklist[ii].desired - tasklist[ii].current;
      double const d_damp(dx.norm() / tasklist[ii].b_max);
      ////dx /= d_damp;		// this is never explicitly stated in [Baerlocher, 2001]
      Vector dx_comp = dx - tasklist[ii].Jacobian * delta_q;
      
      Matrix Jtilda = tasklist[ii].Jacobian * projector;
      
      Matrix Jtilda_inv, delta_projector;
      pseudo_inverse_baerlocher (Jtilda, d_damp, Jtilda_inv, delta_projector,
				 dbgU, dbgsigma, dbgV, dbgdamping);
      
      delta_q = delta_q + Jtilda_inv * dx_comp;
      
      projector = projector - delta_projector;
      
      if (dbgos) {
	*dbgos << dbgpre << "task[" << ii << "]:\n";
	string pre (dbgpre);
	pre += "  ";
	print (dx,                    *dbgos, "dx",              pre);
	print (tasklist[ii].Jacobian, *dbgos, "J",               pre);
	print (dx_comp,               *dbgos, "dx_comp",         pre);
	print (Jtilda,                *dbgos, "Jtilda",          pre);
	print (*dbgU,                 *dbgos, "U",               pre);
	print (*dbgsigma,             *dbgos, "sigma",           pre);
	print (*dbgdamping,           *dbgos, "damping",         pre);
	print (*dbgV,                 *dbgos, "V",               pre);
	print (Jtilda_inv,            *dbgos, "Jtilda_inv",      pre);
	print (delta_projector,       *dbgos, "delta_projector", pre);
	print (delta_q,               *dbgos, "delta_q",         pre);
	print (projector,             *dbgos, "projector",       pre);
	
	Eigen::JacobiSVD<Matrix> svd(projector, Eigen::ComputeFullU | Eigen::ComputeFullV);
	print (svd.singularValues(), *dbgos, "proj. sigma", pre);
	print (svd.matrixU(), *dbgos, "proj. U", pre);
	print (svd.matrixV(), *dbgos, "proj. V", pre);
      }
    }
    
    // // add joint-space "criterion minimization term" (i.e. posture)
    
    // double const xi(-1.0e-4);	// must be negative
    // Vector const & posture_gradient(system.state); // Vector::Ones(ndof)); // some bias...
    
    // delta_q = delta_q + projector * xi * posture_gradient;
    
    if (dbgos) {
      delete dbgU;
      delete dbgV;
      delete dbgsigma;
      delete dbgdamping;
    }
    
    return delta_q;
  }
  
}
