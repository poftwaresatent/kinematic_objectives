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
  
  
  Vector algorithm_unconstrained(vector<Objective *> const & objectives,
				 ostream * dbgos,
				 char const * dbgpre)
  {
    Vector const & xa(objectives[0]->delta_);
    Matrix const & Ja(objectives[0]->Jacobian_);
    Matrix Ja_inv;
    pseudo_inverse_moore_penrose(Ja, Ja_inv);
    Vector ddq(Ja_inv * xa);
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print(xa, *dbgos, "xa", pre);
      print(Ja, *dbgos, "Ja", pre);
      print(Ja_inv, *dbgos, "Ja_inv", pre);
      print(ddq, *dbgos, "ddq", pre);
    }
    
    if (1 == objectives.size()) {
      if (dbgos) {
	*dbgos << dbgpre << "no secondary task\n";
      }
      return ddq;
    }
    
    size_t const ndof(Ja.cols());
    Matrix const Na(Matrix::Identity(ndof, ndof) - Ja_inv * Ja);
    TaskData beta;
    beta.stack(objectives.begin() + 1, objectives.end());
    Vector const & xb(beta.delta_);
    Matrix const & Jb(beta.Jacobian_);
    Matrix Jb_inv;
    pseudo_inverse_moore_penrose(Jb, Jb_inv);
    
    ddq += Na * Jb_inv * xb;
    
    if (dbgos) {
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print(xb, *dbgos, "xb", pre);
      print(Jb, *dbgos, "Jb", pre);
      print(Jb_inv, *dbgos, "Jb_inv", pre);
      print(Na, *dbgos, "Na", pre);
      Matrix tmp;
      tmp = Jb_inv * xb;
      print(tmp, *dbgos, "Jb_inv * xb", pre);
      tmp = Na * Jb_inv;
      print(tmp, *dbgos, "Na * Jb_inv", pre);
      tmp = Na * Jb_inv * xb;
      print(tmp, *dbgos, "Na * Jb_inv * xb", pre);
      print(ddq, *dbgos, "updated ddq", pre);
    }
    
    return ddq;
  }
  
  
  Vector algorithm_constrained(double timestep,
			       Vector const & velocity,
			       Vector const & ddq_unconstrained,
			       vector<Constraint *> const & constraints,
			       ostream * dbgos,
			       char const * dbgpre)
  {
    Vector const & xa(constraints[0]->delta_);
    Matrix const & Ja(constraints[0]->Jacobian_);
    Matrix Ja_inv;
    pseudo_inverse_moore_penrose(Ja, Ja_inv);
    double const t2(timestep * timestep);
    Vector const dq(timestep * velocity);
    Vector ddq((Ja_inv * xa - dq) / t2);
    size_t const ndof(Ja.cols());
    Matrix const Na(Matrix::Identity(ndof, ndof) - Ja_inv * Ja);
    
    if (dbgos) {
      *dbgos << dbgpre << "primary constraint:\n";
      string pre (dbgpre);
      pre += "  ";
      print(xa, *dbgos, "xa", pre);
      print(Ja, *dbgos, "Ja", pre);
      print(Ja_inv, *dbgos, "Ja_inv", pre);
      print(Vector(Ja_inv * xa), *dbgos, "Ja_inv * xa", pre);
      print(dq, *dbgos, "timestep * velocity", pre);
      print(ddq, *dbgos, "ddq", pre);
      print(Na, *dbgos, "Na", pre);
    }
    
    if (1 == constraints.size()) {
      ddq += Na * ddq_unconstrained;
      if (dbgos) {
	*dbgos << dbgpre << "no secondary constraint\n";
	string pre (dbgpre);
	pre += "  ";
	print(Vector(Na * ddq_unconstrained), *dbgos, "Na * ddq_unconstrained", pre);
	print(ddq, *dbgos, "ddq", pre);
      }
      return ddq;
    }
    
    TaskData beta;
    beta.stack(constraints.begin() + 1, constraints.end());
    Vector const & xb(beta.delta_);
    Matrix const & Jb(beta.Jacobian_);
    Matrix Jb_inv;
    pseudo_inverse_moore_penrose(Jb, Jb_inv);
    Matrix const Nb(Matrix::Identity(ndof, ndof) - Jb_inv * Jb);
    
    ddq += Na * ((Jb_inv * xb - dq) / t2 + Nb * ddq_unconstrained);
    
    if (dbgos) {
      *dbgos << dbgpre << "secondary constraint:\n";
      string pre (dbgpre);
      pre += "  ";
      print(xb, *dbgos, "xb", pre);
      print(Jb, *dbgos, "Jb", pre);
      print(Jb_inv, *dbgos, "Jb_inv", pre);
      print(Nb, *dbgos, "Nb", pre);
      print(Vector(Jb_inv * xb), *dbgos, "Jb_inv * xb", pre);
      print(Vector((Jb_inv * xb - dq) / t2), *dbgos, "(Jb_inv * xb - dq) / t2", pre);
      print(Vector(Nb * ddq_unconstrained), *dbgos, "Nb * ddq_unconstrained", pre);
      print(Vector((Jb_inv * xb - dq) / t2 + Nb * ddq_unconstrained), *dbgos, "(Jb_inv * xb - dq) / t2 + Nb * ddq_unconstrained", pre);
      print(Vector(Na * ((Jb_inv * xb - dq) / t2 + Nb * ddq_unconstrained)), *dbgos, "Na * ((Jb_inv * xb - dq) / t2 + Nb * ddq_unconstrained)", pre);
      print(ddq, *dbgos, "updated ddq", pre);
    }
    
    return ddq;
  }
  
}
