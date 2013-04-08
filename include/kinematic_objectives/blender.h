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
#include <list>


namespace kinematic_objectives {
  
  class CompoundObjective;
  
  /**
     The algorithm which takes a compound objective and turns it into
     an overall update. Slated to become an abstract base class at
     some point, and provide generic inspection (introspection?)
     capabilities so that users get informaed why some aspect of the
     motion they expect is not getting fulfilled (e.g. "why is the
     tray tilting although I had told the robot to keep it horizontal"
     --> "it is because of an interference between a joint limit and
     obstacle avoidance").
     
     The Blender should provide some quantitative "meta feedback"
     about its results. For example, which unilateral constraints were
     activated, the resulting constraint nullspace, whether that
     created any conflicts with hard objectives, whether there were
     any algorithmic singularities between constraints and hard
     objectives. For soft objectives, it should at least be possible
     to specify the remaining space of motion freedom they had
     available collectively.
     
     \todo [medium] make the attributes protected or private; [medium]
     extract a nice base class and provide at least two implementation
     (based on Chiaverini:1997 and Baerlocher:2001, maybe also one for
     Siciliano:1991); [low] find a better name.
  */  
  class Blender
  {
  public:
    typedef list<CompoundObjective *> path_t;
    
    Blender(double timestep, ostream * dbgos, string const & dbgpre);
    
    virtual ~Blender();
    
    virtual void clear();
    
    virtual void update();
    
    /**
       \note This is still under somewhat active development, see
       comments in the implementation. Inspired by a combination of
       Chiaverini:1997 and Baerlocher:2001 with the aim of achieving
       proper decoupling between priority levels while also switching
       unilateral constraints (e.g. joint limits and obstacle
       avoidance) on and off based on the current joint positions and
       velocities. The latter is an important feature, but the
       implementation is rather more convoluted than hoped.
    */
    virtual void updateCompoundObjective(CompoundObjective * wpt);
    
    ostream * dbgos_;
    string dbgpre_;
    string dbgpre2_;
    
    double timestep_;
    path_t path_;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_BLENDER_HPP
