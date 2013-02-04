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

#ifndef KINEMATIC_ELASTIC_JOINT_LIMIT_CONSTRAINT_HPP
#define KINEMATIC_ELASTIC_JOINT_LIMIT_CONSTRAINT_HPP

#include "task.hpp"


namespace kinematic_elastic {
  
  class JointLimitConstraint
    : public Task
  {
  public:
    void init(size_t ndof);
    
    /**
       For every joint that lies outside hard joint limits, create a
       task dimension that brings it back to the corresponding soft
       limit. All such dimensions are stacked into our TaskData
       fields. Furthermore, the list of violating joint indices is
       placed in locked_joints_. If there are no joints that violate
       their hard joint limits, the task and locked list will be
       empty.
    */
    virtual void update(Model const & model);
    
    virtual bool isActive() const;
    
    /**
       Nx4 matrix, one row per joint, where
       - col[0] is the lower hard limit
       - col[1] is the lower soft limit
       - col[2] is the upper soft limit
       - col[3] is the upper hard limit
    */
    Matrix limits_;
    
    vector<size_t> locked_joints_;
  };
  
}

#endif // KINEMATIC_ELASTIC_JOINT_LIMIT_CONSTRAINT_HPP
