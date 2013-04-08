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

#include "interactive_blender.h"


namespace kinematic_objectives {

  namespace demo {
    
    
    InteractiveBlender::
    InteractiveBlender(Blender * blender,
		       ostream * dbgos,
		       string const & dbgpre)
      : ee_      (0.2, 0.0, 1.0, 0.0, 0.5),
	ee_ori_  (0.1, 0.0, 1.0, 0.0, 0.3),
	z_angle_ (0),
	base_    (0.2, 0.0, 1.0, 0.5, 0.5),
	repulsor_(1.5, 1.0, 0.5, 0.0, 0.2),
	obstacle_(1.5, 0.7, 0.0, 0.2, 0.5),
	robot_(*this,
	       repulsor_,
	       z_angle_,
	       &(ee_.point_),
	       &(base_.point_)),
	blender_(blender)
    {
      handles_.push_back(&ee_);
      handles_.push_back(&ee_ori_);
      handles_.push_back(&base_);
      handles_.push_back(&repulsor_);
      handles_.push_back(&obstacle_);
    }
    
    
    void InteractiveBlender::
    init(Vector const & state)
    {
      // vector<NormalCompoundObjective *> wpt;
      // for (size_t ii(0); ii < 10; ++ii) {
      // 	wpt.push_back(new NormalCompoundObjective(*this, repulsor_, z_angle_));
      // }
      // wpt[0]->setNeighbors(start, wpt[1]);
      // for (size_t ii(1); ii < wpt.size() - 1; ++ii) {
      // 	wpt[ii]->setNeighbors(wpt[ii-1], wpt[ii+1]);
      // }
      // wpt[wpt.size() - 1]->setNeighbors(wpt[wpt.size() - 2], goal);
      
      robot_.init(state, Vector::Zero(state.size()));
    }
      
      
    void InteractiveBlender::
    update()
    {
      z_angle_ = atan2(ee_ori_.point_[1] - ee_.point_[1], ee_ori_.point_[0] - ee_.point_[0]);
      blender_->update(&robot_);
    }
    
    
    void InteractiveBlender::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      robot_.draw(cr, weight, pixelsize);
      
      for (size_t ii(0); ii < handles_.size(); ++ii) {
	handles_[ii]->draw(cr, weight, pixelsize);
      }
      
      cairo_set_source_rgba(cr, 0.0, 1.0, 0.0, 0.5);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr, ee_.point_[0], ee_.point_[1]);
      cairo_line_to(cr, ee_ori_.point_[0], ee_ori_.point_[1]);
      cairo_stroke(cr);
    }
    
  }
  
}
