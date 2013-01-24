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
#include <Eigen/SVD>
#include <iostream>
#include <deque>


namespace kinematic_elastic {
  
  
  Vector algorithm (Model const & model,
		    Vector const & state,
		    tasklist_t const & tasklist,
		    ostream * dbgos,
		    char const * dbgpre)
  {
    deque<task_s const *> tl(tasklist.size());
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      tl[ii] = &tasklist[ii];
    }
    
    task_s jl;
    model.createJointLimitTask(state, jl);
    if (jl.ndim) {
      tl.push_front(&jl);
    }
    
    size_t const ndof(state.size());
    Vector dq(Vector::Zero(ndof));
    
    if (tl.empty()) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (no tasks):\n";
	print (dq, *dbgos, "dq", string(dbgpre) + "  ");
      }
      return dq;
    }
    
    Vector dxa(tl[0]->desired - tl[0]->current);
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (tl[0]->Jacobian, Ja_inv);
    dq += Ja_inv * dxa;
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (dxa,             *dbgos, "dxa",    pre);
      print (tl[0]->Jacobian, *dbgos, "Ja",     pre);
      print (Ja_inv,          *dbgos, "Ja_inv", pre);
      print (dq,              *dbgos, "dq",     pre);
    }
    
    if (tl.size() == 1) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (single task):\n";
      }
      return dq;
    }
    
    Vector dxb(tl[0]->ndim + tl[1]->ndim);
    dxb.block(0,           0, tl[0]->ndim, 1) = dxa;
    dxb.block(tl[0]->ndim, 0, tl[1]->ndim, 1) = tl[1]->desired - tl[1]->current;
    Matrix Jb(tl[0]->ndim + tl[1]->ndim, ndof);
    Jb.block(0,           0, tl[0]->ndim, ndof) = tl[0]->Jacobian;
    Jb.block(tl[0]->ndim, 0, tl[1]->ndim, ndof) = tl[1]->Jacobian;
    Matrix Jb_inv;
    double tmp(1.0 / (tl[0]->b_max < tl[1]->b_max ? tl[0]->b_max : tl[1]->b_max));
    pseudo_inverse_damped (Jb, tmp, Jb_inv);
    Matrix Na(Matrix::Identity(ndof, ndof) - Ja_inv * tl[0]->Jacobian);
    dq +=  Na * Jb_inv * dxb;

    if (dbgos) {
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Na,     *dbgos, "Na",     pre);
      print (dxb,    *dbgos, "dxb",    pre);
      print (Jb,     *dbgos, "Jb",     pre);
      print (Jb_inv, *dbgos, "Jb_inv", pre);
      print (dq,     *dbgos, "dq",     pre);
    }
    
    if (tl.size() == 2) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (two tasks):\n";
      }
      return dq;
    }
    
    Matrix Nb(Matrix::Identity(ndof, ndof) - Jb_inv * Jb);
    size_t nrest(0);
    for (size_t it(2); it < tl.size(); ++it) {
      nrest += tl[it]->ndim;
    }
    Vector dxc(nrest);
    Matrix Jc(nrest, ndof);
    tmp = numeric_limits<double>::max();
    for (size_t it(2), ir(0); it < tl.size(); ir += tl[it++]->ndim) { // obscure update expression, grr
      dxc.block(ir, 0, tl[it]->ndim, 1) = tl[it]->desired - tl[it]->current;
      Jc.block( ir, 0, tl[it]->ndim, ndof) = tl[it]->Jacobian;
      if (tl[it]->b_max < tmp) {
	tmp = tl[it]->b_max; // beware: tmp is 1.0/lambda here (whereas it was lambda above)
      }
    }
    Matrix Jc_inv;
    pseudo_inverse_damped (Jc, 1.0 / tmp, Jc_inv);
    dq += Na * Nb * Jc_inv * dxc;
    
    if (dbgos) {
      *dbgos << dbgpre << "remainder:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Nb,     *dbgos, "Nb",     pre);
      print (dxc,    *dbgos, "dxc",    pre);
      print (Jc,     *dbgos, "Jc",     pre);
      print (Jc_inv, *dbgos, "Jc_inv", pre);
      print (dq,     *dbgos, "dq",     pre);
      *dbgos << dbgpre << "DONE\n";
    }
    
    return dq;
  }
  
}
