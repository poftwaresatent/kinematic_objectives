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

#ifndef KINEMATIC_OBJECTIVES_CONSTRAINT_TELEPORTING_BLENDER_HPP
#define KINEMATIC_OBJECTIVES_CONSTRAINT_TELEPORTING_BLENDER_HPP

#include <kinematic_objectives/blender.h>


namespace kinematic_objectives {
  
  /**
     An experimental (but working) implementation of the Blender
     interface.
  */  
  class ConstraintTeleportingBlender
    : public Blender
  {
  public:
    ConstraintTeleportingBlender(double timestep, ostream * dbgos, string const & dbgpre);
    
    /**
       \note This is still somewhat experimental. See
       comments in the implementation. Inspired by a combination of
       Chiaverini:1997 and Baerlocher:2001 with the aim of achieving
       proper decoupling between priority levels while also switching
       unilateral constraints (e.g. joint limits and obstacle
       avoidance) on and off based on the current joint positions and
       velocities. The latter is an important feature, but the
       implementation is rather more convoluted than hoped.
    */
    virtual void update(CompoundObjective * wpt);
    
    ostream * dbgos_;
    string dbgpre_;
    string dbgpre2_;
    
    double timestep_;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_CONSTRAINT_TELEPORTING_BLENDER_HPP
