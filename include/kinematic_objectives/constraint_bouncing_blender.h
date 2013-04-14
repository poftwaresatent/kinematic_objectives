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

#ifndef KINEMATIC_OBJECTIVES_CONSTRAINT_BOUNCING_BLENDER_HPP
#define KINEMATIC_OBJECTIVES_CONSTRAINT_BOUNCING_BLENDER_HPP

#include <kinematic_objectives/blender.h>


namespace kinematic_objectives {
  
  /**
     A blender based on the "classical" approach [Siciliano:1991]
     which handles unilateral constraints in a straightforward manner:
     it computes their desired displacement, then updates the
     kinematic model by weighting that displacement and projecting the
     current velocities into the constraint nullspace. The effect is
     that unlateral constraints do get switched on, by they typically
     stay that way because the constraint nullspace keeps the
     objectives from pulling the state away from the violation. Also,
     it has a tendency to bounce off of constraints, due to (as far as
     I can tell at this moment) a combination of linearization and
     discretization errors inherent in the approach.  Thus it should
     be clear that this blender is useful mostly for development and
     testing.
  */
  class ConstraintBouncingBlender
    : public Blender
  {
  public:
    explicit ConstraintBouncingBlender(double stepsize, ostream * dbgos, string const & dbgpre);
    
    virtual void update(KinematicModel & model, CompoundObjective * wpt);
  };
  
}

#endif // KINEMATIC_OBJECTIVES_CONSTRAINT_BOUNCING_BLENDER_HPP
