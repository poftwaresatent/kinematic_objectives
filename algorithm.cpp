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
  
  
  /**
     \pre alpha must be non-null and non-singular.
  */
  static Vector compute_dq (TaskData * alpha,
			    TaskData * beta,
			    TaskData * gamma,
			    ostream * dbgos,
			    char const * dbgpre)
  {
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (alpha->Jacobian_, Ja_inv);
    Vector dq(Ja_inv * alpha->delta_);
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (alpha->delta_, *dbgos, "dxa",    pre);
      print (alpha->Jacobian_, *dbgos, "Ja",     pre);
      print (Ja_inv,    *dbgos, "Ja_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
    }
    
    if (0 == beta) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (single task)\n";
      }
      return dq;
    }
    
    Matrix Jb_inv;
    pseudo_inverse_damped (beta->Jacobian_, 1.0 / beta->step_hint_, Jb_inv);
    size_t const ndof(alpha->Jacobian_.cols());
    Matrix Na(Matrix::Identity(ndof, ndof) - Ja_inv * alpha->Jacobian_);
    dq +=  Na * Jb_inv * beta->delta_;
    
    if (dbgos) {
      Matrix Na_times_Jb_inv(Na * Jb_inv);
      Vector dqb(Na_times_Jb_inv * beta->delta_);
      
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Na,              *dbgos, "Na", pre);
      print (beta->delta_,       *dbgos, "dxb", pre);
      print (beta->Jacobian_,       *dbgos, "Jb", pre);
      print (Jb_inv,          *dbgos, "Jb_inv", pre);
      print (Na_times_Jb_inv, *dbgos, "Na * Jb_inv", pre);
      print (dqb,             *dbgos, "Na * Jb_inv * dxb", pre);
      print (dq,              *dbgos, "dq", pre);
    }
    
    if (0 == gamma) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (two tasks)\n";
      }
      return dq;
    }
    
    Matrix Nb(Matrix::Identity(ndof, ndof) - Jb_inv * beta->Jacobian_);
    Matrix Jc_inv;
    pseudo_inverse_damped (gamma->Jacobian_, 1.0 / gamma->step_hint_, Jc_inv);
    dq += Na * Nb * Jc_inv * gamma->delta_;
    
    if (dbgos) {
      *dbgos << dbgpre << "remainder:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Nb,        *dbgos, "Nb",     pre);
      print (gamma->delta_, *dbgos, "dxc",    pre);
      print (gamma->Jacobian_, *dbgos, "Jc",     pre);
      print (Jc_inv,    *dbgos, "Jc_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
      *dbgos << dbgpre << "DONE\n";
    }
    
    return dq;
  }
  
}
