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

#ifndef KINEMATIC_OBJECTIVES_WAYPOINT_HPP
#define KINEMATIC_OBJECTIVES_WAYPOINT_HPP

#include <kinematic_objectives/types.h>


namespace kinematic_objectives {
  
  class KinematicModel;
  class Objective;
  
  
  /**
     A collection of objectives that can be turned into a prioritized
     and/or summed overall update using a blender. Currently hardcodes
     the notion of "constraint" (delta-position-based update with
     corresponding nullspace), "hard objective" (a.k.a. task in the
     controls community, here implemented as an acceleration-based
     update with corresponding nullspace), and "soft objective" (which
     is acceleration based and provides the somewhat novel way of
     pretending we have virtual forces in the absence of system
     dynamic models).
     
     \todo [high] generalize this (away from the old "waypoint" idea
     which is still present in the demo code); [medium] make the
     attributes protected or private.
   */  
  class CompoundObjective
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    virtual ~CompoundObjective();
    
    virtual void init(KinematicModel const & model);
    
    vector<Objective *> unilateral_constraints_;
    vector<Objective *> hard_objectives_;
    vector<Objective *> soft_objectives_;
  };  
  
}

#endif // KINEMATIC_OBJECTIVES_WAYPOINT_HPP
