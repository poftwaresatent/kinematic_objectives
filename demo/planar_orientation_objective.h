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

#ifndef KINEMATIC_OBJECTIVES_DEMO_PLANAR_ORIENTATION_ORIENTATION_OBJECTIVE_HPP
#define KINEMATIC_OBJECTIVES_DEMO_PLANAR_ORIENTATION_ORIENTATION_OBJECTIVE_HPP

#include <kinematic_objectives/objective.h>

namespace kinematic_objectives {
  
  namespace demo {
    
    /**
       \todo [low] attributes should be protected or private
    */
    class PlanarOrientationObjective
      : public Objective
    {
    public:
      PlanarOrientationObjective(string const & name,
				 size_t node,
				 double kp,
				 double kd);
      
      virtual void init(KinematicModel const & model);
      virtual void update(KinematicModel const & model);
      virtual double computeResidualErrorMagnitude(Vector const & ee) const;
      
      double kp_;
      double kd_;
      size_t node_;
      double angle_;
      double goal_;
    };
    
  }
  
}

#endif // KINEMATIC_OBJECTIVES_DEMO_PLANAR_ORIENTATION_ORIENTATION_OBJECTIVE_HPP
