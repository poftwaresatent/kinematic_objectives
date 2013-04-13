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

#ifndef KINEMATIC_OBJECTIVES_ACHIEVABILITY_HPP
#define KINEMATIC_OBJECTIVES_ACHIEVABILITY_HPP

#include <kinematic_objectives/types.h>


namespace kinematic_objectives {
  
  class KinematicModel;
  class Objective;
  class CompoundObjective;
  
  class Achievability
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    typedef enum {
      UNILATERAL_CONSTRAINT,
      HARD_OBJECTIVE,
      SOFT_OBJECTIVE
    } objective_tag_t;
    
    objective_tag_t tag_;
    Objective const * objective_;
    size_t available_dimension_;
    size_t required_dimension_;
    Vector residual_error_;
    double residual_error_magnitude_;
    
    Vector nullspace_residuals_;
    
    // Matrix jbar_inverse_;
    // size_t jbar_range_;
    // requires SVD of J: Matrix nullspace_residuals_;
    // requires SVD of J: Matrix removed_directions_;
    
    /**
       Needs an extra SVD for each active objective. So don't call
       this too often... or if you don't care about the nullspace
       residuals, make a lightweight copy of this function which skips
       that step.
       
       \note Soft objectives are a bit special because they do not
       individually get projected. So here we compute a few spurious
       SVDs just in order to find out how many dimensions each soft
       objective wants and gets. Assuming that this info is not
       computed too often, that should not result in an appreciable
       slowdown.
    */

    static void compute(KinematicModel & model,
			CompoundObjective const & co,
			vector<Achievability> & information);
    
    static void print(vector<Achievability> const & information,
		      ostream & os, string const & pfx);
  };
  
}

#endif // KINEMATIC_OBJECTIVES_ACHIEVABILITY_HPP
