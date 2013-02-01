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
			   tasklist_t const & tl,
			   ostream * dbgos,
			   char const * dbgpre)
  {
    Vector dxb;
    stack(tl[0]->dx, tl[1]->dx, dxb);
    
    Matrix Jb;
    stack(tl[0]->Jx, tl[1]->Jx, Jb);
    
    // could add support for more than 2 tasks later...
    
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (tl[0]->Jx, Ja_inv);
    
    size_t const ndof(state.size());
    Matrix Na = Matrix::Identity(ndof, ndof) - Ja_inv * tl[0]->Jx;
    
    Matrix Jb_inv;
    double foo(tl[0]->b_max); // well... how to find a good lambda?
    if (tl[1]->b_max < foo) {
      foo = tl[1]->b_max;
    }
    pseudo_inverse_damped (Jb, 1.0 / foo, Jb_inv);
    
    // would need Nb for supporting more than 2 tasks... later...
    
    Vector dq = Ja_inv * tl[0]->dx + Na * Jb_inv * dxb;
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print (tl[0]->dx, *dbgos, "dxa",    pre);
      print (tl[0]->Jx, *dbgos, "Ja",     pre);
      print (Ja_inv,    *dbgos, "Ja_inv", pre);
      print (Na,        *dbgos, "Na",     pre);
      print (dxb,       *dbgos, "dxb",    pre);
      print (Jb,        *dbgos, "Jb",     pre);
      print (Jb_inv,    *dbgos, "Jb_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
    }
    
    return dq;
  }
  
}
