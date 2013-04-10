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

#ifndef KINEMATIC_OBJECTIVES_BLENDER_HPP
#define KINEMATIC_OBJECTIVES_BLENDER_HPP

#include <kinematic_objectives/types.h>


namespace kinematic_objectives {
  
  class CompoundObjective;
  
  /**
     Interface for an algorithm that takes a compound objective and
     turns it into an overall update. Provides generic inspection
     (introspection?)  capabilities to inform users when and why some
     aspect of the motion they expect is not getting fulfilled
     (e.g. "the tray is tilting because it gets overridden by an
     interference between joint limits and obstacle avoidance"). The
     actual feedback data structure lives in CompoundObjective (at the
     time of writing).
     
     \todo [high] create a blender subclass based on Chiaverini:1997;
     [low] make the attributes protected or private; [low] maybe find
     a different name.
  */  
  class Blender
  {
  public:
    string const name_;
    
    explicit Blender(string const & name)
      : name_(name) { }
    
    virtual ~Blender()
    { /* nop */ }
    
    virtual void update(CompoundObjective * wpt) = 0;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_BLENDER_HPP
