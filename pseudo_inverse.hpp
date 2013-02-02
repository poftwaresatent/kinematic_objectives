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

#ifndef KINEMATIC_ELASTIC_PSEUDO_INVERSE_HPP
#define KINEMATIC_ELASTIC_PSEUDO_INVERSE_HPP

#include "kinematic_elastic.hpp"

namespace kinematic_elastic {
  
  void pseudo_inverse_nonsingular(Matrix const & mx,
				  Matrix & inv);
  
  void pseudo_inverse_moore_penrose(Matrix const & mx,
				    Matrix & inv);
  
  void pseudo_inverse_damped(Matrix const & mx,
			     double lambda,
			     Matrix & inv);
  
  /**
     Uses the scheme in [Baerlocher 2001] to determine damping factor
     lambda based on the smallest sigma and the given d_damp
     parameter. See equations (4.17) and figure 4.7. He computes
     d_damp as dx.norm()/b_max where b_max is a fixed parameter.
  */
  void pseudo_inverse_baerlocher(Matrix const & mx,
				 double d_damp,
				 Matrix & inv,
				 Matrix & delta_projector,
				 Matrix * dbgU,
				 Vector * dbgsigma,
				 Matrix * dbgV,
				 Vector * dbgdamping);
  
}

#endif // KINEMATIC_ELASTIC_PSEUDO_INVERSE_HPP
