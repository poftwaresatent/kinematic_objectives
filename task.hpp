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

#ifndef KINEMATIC_ELASTIC_TASK_HPP
#define KINEMATIC_ELASTIC_TASK_HPP

#include "kinematic_elastic.hpp"
#include <limits>


namespace kinematic_elastic {
  
  class Model;
  
  
  class TaskData
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    TaskData();
    
    void stack(TaskData const & t1, TaskData const & t2);
    
    template<typename iterator_t>
    void stack(iterator_t begin, iterator_t end)
    {
      size_t ttnrows(0);
      step_hint_ = numeric_limits<double>::max();
      for (iterator_t ii(begin); ii != end; ++ii) {
	ttnrows += (*ii)->Jacobian_.rows();
	if (step_hint_ > (*ii)->step_hint_) {
	  step_hint_ = (*ii)->step_hint_;
	}
      }
      size_t const ndof((*begin)->Jacobian_.cols());
      delta_.resize(ttnrows);
      Jacobian_.resize(ttnrows, ndof);
      for (size_t row(0); begin != end; row += (*begin++)->Jacobian_.rows()) {
	delta_.block(   row, 0, (*begin)->Jacobian_.rows(),    1) = (*begin)->delta_;
	Jacobian_.block(row, 0, (*begin)->Jacobian_.rows(), ndof) = (*begin)->Jacobian_;
      }
    }
    
    Vector delta_;		// desired - current
    Matrix Jacobian_;
    double step_hint_;
  };
  
  
  class Task
    : public TaskData
  {
  public:
    virtual ~Task();
    
    virtual bool init(Model const & model) = 0;
    virtual bool update(Model const & model) = 0;
  };
  
  
  class Objective
    : public Task
  {
  public:
    Objective();
    
    void configure(double kp, double kd); // XXXX to do: add maxvel and saturation method
    
    Vector computePD(Vector const & pos_act,
		     Vector const & vel_act,
		     Vector const & pos_des) const; // XXXX to do: maybe allow vel_des also?
    
    double kp_;
    double kd_;
  };
  
}

#endif // KINEMATIC_ELASTIC_TASK_HPP
