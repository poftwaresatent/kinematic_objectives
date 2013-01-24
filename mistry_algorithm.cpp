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

#include "mistry_algorithm.hpp"
#include "pseudo_inverse.hpp"
#include "print.hpp"
#include <Eigen/SVD>
#include <iostream>


namespace kinematic_elastic {
  
  
  Vector mistry_algorithm (Model const & model,
			   Vector const & state,
			   tasklist_t const & tasklist,
			   ostream * dbgos,
			   char const * dbgpre)
  {
    Vector const dxa(tasklist[0].desired - tasklist[0].current);
    
    Matrix const & Ja(tasklist[0].Jacobian);
    
    Vector dxb(tasklist[0].ndim + tasklist[1].ndim);
    dxb.block(0,                0, tasklist[0].ndim, 1) = dxa;
    dxb.block(tasklist[0].ndim, 0, tasklist[1].ndim, 1) = tasklist[1].desired - tasklist[1].current;
    
    size_t const ndof(state.size());
    Matrix Jb(tasklist[0].ndim + tasklist[1].ndim, ndof);
    Jb.block(0,                0, tasklist[0].ndim, ndof) = Ja;
    Jb.block(tasklist[0].ndim, 0, tasklist[1].ndim, ndof) = tasklist[1].Jacobian;
    
    // could add support for more than 2 tasks later...
    
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (Ja, Ja_inv);
    
    Matrix Na = Matrix::Identity(ndof, ndof) - Ja_inv * Ja;
    
    Matrix Jb_inv;
    double foo(tasklist[0].b_max); // well... how to find a good lambda?
    if (tasklist[1].b_max < foo) {
      foo = tasklist[1].b_max;
    }
    pseudo_inverse_damped (Jb, 1.0 / foo, Jb_inv);
    
    // would need Nb for supporting more than 2 tasks... later...
    
    Vector dq = Ja_inv * dxa + Na * Jb_inv * dxb;
    return dq;
  }
  
}
