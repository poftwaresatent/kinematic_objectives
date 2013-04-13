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

#ifndef KINEMATIC_OBJECTIVES_PSEUDO_INVERSE_HPP
#define KINEMATIC_OBJECTIVES_PSEUDO_INVERSE_HPP

#include <kinematic_objectives/types.h>

namespace kinematic_objectives {
  
  /**
     For well-behaved matrices, \f$ M^{\#} = M^T (M M^T) ^{-1} \f$.
  */
  void pseudo_inverse_nonsingular(Matrix const & mx,
				  Matrix & inv);
  
  /**
     In case you need to subsequently inspect what went on during
     matrix inversion.
  */
  struct PseudoInverseFeedback {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Vector original_sigma;
    Vector regularized_sigma;
    Matrix output_space;	/**< a.k.a. "matrix U" */
    Matrix input_space;		/**< a.k.a. "matrix V" */
  };
  
  /**
     \note I think this ended up being the only pseudo-inverse
     implementation actually used by the blender. I'm keeing the
     others around though for the (hopefully soon) split of Blender
     into a generic superclass with a handfull of specific
     implementations based on existing state of the art and novel
     developments.
     
     \todo [high] contains a magic number (prec = 1e-3) that gets
     converted into an eigenvaluee threshold.
  */
  void pseudo_inverse_moore_penrose(/**
				       Matrix to be inverted
				    */
				    Matrix const & mx,
				    /**
				       Inverse computed via thresholded SVD.
				    */
				    Matrix & inv,
				    /**
				       Optional output parameter:
				       nullspace projection update
				       (N_proj -= *dproj)
				    */
				    Matrix * dproj = 0,
				    /**
				       Optional output parameter.
				    */
				    PseudoInverseFeedback * fb = 0);
  
  void pseudo_inverse_damped(Matrix const & mx,
			     double lambda,
			     Matrix & inv);
  
  /**
     Uses the scheme in [Baerlocher:2001] to determine damping factor
     lambda based on the smallest sigma and the given d_damp
     parameter. See equations (4.17) and figure 4.7. He computes
     d_damp as dx.norm()/b_max where b_max is a fixed parameter.
  */
  void pseudo_inverse_baerlocher(Matrix const & mx,
				 double d_damp,
				 Matrix & inv,
				 Matrix & bias_projector,
				 Matrix * dbgU,
				 Vector * dbgsigma,
				 Matrix * dbgV,
				 Vector * dbgdamping);
  
}

#endif // KINEMATIC_OBJECTIVES_PSEUDO_INVERSE_HPP
