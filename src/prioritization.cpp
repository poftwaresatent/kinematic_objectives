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

#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  void prioritization_siciliano1991(Matrix const & N_init,
				    vector<Objective*> const & objectives,
				    Vector & bias_res,
				    Matrix & N_res,
				    ostream * dbgos,
				    string const & dbgpre)
  {
    bias_res = Vector::Zero(N_init.rows());
    N_res = N_init;
    
    Matrix Jbinv;		// pseudo-inverse of J_bar (which is J * N)
    Matrix Nup;			// nullspace updater: N -= N_up at each hierarchy level
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "initial nullspace", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of nullspace", pre);
    }
    
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      
      if ( ! objectives[ii]->isActive()) {
	if (dbgos) {
	  *dbgos << dbgpre << "objective " << ii << " is INACTIVE\n";
	}	
	continue;
      }
      
      if (dbgos) {
	*dbgos << dbgpre << "objective " << ii << " \"" << objectives[ii]->name_ << "\":\n";
	string pre (dbgpre);
	pre += "  ";
	print(objectives[ii]->getBias(), *dbgos, "objective bias", pre);
	print(objectives[ii]->getJacobian(), *dbgos, "objective Jacobian", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(objectives[ii]->getJacobian(), Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of objective Jacobian", pre);
	Vector vtmp;
	vtmp = objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res;
	print(vtmp, *dbgos, "bias_comp", pre);
	Matrix mtmp;
	mtmp = objectives[ii]->getJacobian() * N_res;
	print(mtmp, *dbgos, "J_bar", pre);
      }
      
      pseudo_inverse_moore_penrose(objectives[ii]->getJacobian() * N_res, Jbinv, &Nup,
				   &objectives[ii]->jbar_svd_);
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(objectives[ii]->jbar_svd_.singular_values, *dbgos, "eigenvalues of J_bar", pre);
	print(Jbinv, *dbgos, "J_bar pseudo inverse", pre);
	print(Nup, *dbgos, "nullspace projector update", pre);
	Vector vtmp;
        vtmp = Jbinv * (objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res);
	print(vtmp, *dbgos, "bias update", pre);
        vtmp = objectives[ii]->getJacobian() * vtmp;
	print(vtmp, *dbgos, "biased objective update", pre);
        vtmp = objectives[ii]->getJacobian() * Jbinv * objectives[ii]->getBias();
	print(vtmp, *dbgos, "unbiased objective update", pre);
      }
      
      bias_res += Jbinv * (objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res);
      
      if (Nup.cols() == 0) {
	if (dbgos) {
	  *dbgos << dbgpre << "NO DEGREES OF FREEDOM LEFT\n";
	}
	break;
      }
      
      N_res -= Nup;
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(bias_res, *dbgos, "accumulated bias", pre);
	print(N_res, *dbgos, "accumulated nullspace projector", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of nullspace projector", pre);
      }
      
    }
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "final nullspace projectir", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of final nullspace projector", pre);
    }

  }
  
}
