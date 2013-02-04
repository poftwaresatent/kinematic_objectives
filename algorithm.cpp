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

#include "algorithm.hpp"
#include "pseudo_inverse.hpp"
#include "print.hpp"
#include "task.hpp"
#include "joint_limits.hpp"
#include <Eigen/LU>
#include <iostream>


namespace kinematic_elastic {

  void perform_prioritization(Matrix const & N_init,
			      vector<Task*> const & tasks,
			      Vector & delta_res,
			      Matrix & N_res,
			      ostream * dbgos,
			      char const * dbgpre)
  {
    delta_res = Vector::Zero(N_init.rows());
    N_res = N_init;
    
    Matrix Jbinv;
    Matrix Nup;
    for (size_t ii(0); ii < tasks.size(); ++ii) {
      
      if (dbgos) {
	*dbgos << dbgpre << "task " << ii << ":\n";
	string pre (dbgpre);
	pre += "  ";
	print(tasks[ii]->delta_, *dbgos, "task delta", pre);
	print(tasks[ii]->Jacobian_, *dbgos, "task Jacobian", pre);
	Vector vtmp;
	vtmp = tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res;
	print(vtmp, *dbgos, "delta_comp", pre);
	Matrix mtmp;
	mtmp = tasks[ii]->Jacobian_ * N_res;
	print(mtmp, *dbgos, "J_bar", pre);
      }
      
      pseudo_inverse_moore_penrose(tasks[ii]->Jacobian_ * N_res, Jbinv, &Nup);
      delta_res += Jbinv * (tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res);
      N_res -= Nup;
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(Jbinv, *dbgos, "J_bar pseudo inverse", pre);
	print(Nup, *dbgos, "nullspace update", pre);
	Vector vtmp;
        vtmp = Jbinv * (tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res);
	print(vtmp, *dbgos, "delta update", pre);
	print(delta_res, *dbgos, "accumulated delta", pre);
	print(N_res, *dbgos, "accumulated nullspace", pre);
      }
      
    }
  }
  
}
