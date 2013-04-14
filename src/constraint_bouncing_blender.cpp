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

#include <kinematic_objectives/constraint_bouncing_blender.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>
#include <kinematic_objectives/integrator.h>


namespace kinematic_objectives {
  
  
  ConstraintBouncingBlender::
  ConstraintBouncingBlender(Integrator const * integrator)
    : Blender("ConstraintBouncingBlender", integrator)
  {
  }
  
  
  void ConstraintBouncingBlender::
  update(KinematicModel & model, CompoundObjective * wpt)
  {
    for (size_t ii(0); ii < wpt->unilateral_constraints_.size(); ++ii) {
      wpt->unilateral_constraints_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(model);
    }
    
    ssize_t const ndof(model.getJointPosition().size());
    Vector bias_c;
    Matrix np_c;
    prioritization_.projectObjectives(Matrix::Identity(ndof, ndof),
				      Vector::Zero(ndof),
				      wpt->unilateral_constraints_,
				      bias_c,
				      np_c);
    
    Vector bias_h;
    Matrix np_h;
    prioritization_.projectObjectives(np_c,
				      bias_c,
				      wpt->hard_objectives_,
				      bias_h,
				      np_h);
    
    Vector bias_s;
    prioritization_.addUpObjectives(np_h,
				    // ////Vector::Zero(ndof),
				    // bias_h,
				    wpt->soft_objectives_,
				    bias_s);
    
    Vector q_next, qd_next;
    integrator_->compute(bias_c + bias_h + bias_s,
			 model.getJointPosition(),
			 model.getJointVelocity(),
			 q_next,
			 qd_next);
    
    model.update(q_next, qd_next);
  }
  
}
