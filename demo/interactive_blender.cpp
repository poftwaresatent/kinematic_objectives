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
    InteractiveBlender(double timestep,
		       ostream * dbgos,
		       string const & dbgpre)
      : Blender(timestep, dbgos, dbgpre),
	eestart_(0.2, 0.0, 1.0, 0.0, 0.5),
	eestartori_(0.1, 0.0, 1.0, 0.0, 0.3),
	z_angle_(0),
	basestart_(0.2, 0.0, 1.0, 0.5, 0.5),
	eegoal_(0.2, 0.0, 0.0, 1.0, 0.5),
	basegoal_(0.2, 0.0, 0.5, 1.0, 0.5),
	repulsor_(1.5, 1.0, 0.5, 0.0, 0.2),
	obstacle_(1.5, 0.7, 0.0, 0.2, 0.5)
    {
      handles_.push_back(&eestart_);
      handles_.push_back(&eestartori_);
      handles_.push_back(&basestart_);
      handles_.push_back(&eegoal_);
      handles_.push_back(&basegoal_);
      handles_.push_back(&repulsor_);
      handles_.push_back(&obstacle_);
    }
    
    
    void InteractiveBlender::
    init(Vector const & state)
    {
      clear();
      
      BoundaryCompoundObjective * start(new BoundaryCompoundObjective(*this,
								      repulsor_,
								      z_angle_,
								      &(eestart_.point_),
								      &(basestart_.point_)));
      BoundaryCompoundObjective * goal(new BoundaryCompoundObjective(*this,
								     repulsor_,
								     z_angle_,
								     &(eegoal_.point_),
								     &(basegoal_.point_)));
      vector<NormalCompoundObjective *> wpt;
      for (size_t ii(0); ii < 10; ++ii) {
	wpt.push_back(new NormalCompoundObjective(*this, repulsor_, z_angle_));
      }
      wpt[0]->setNeighbors(start, wpt[1]);
      for (size_t ii(1); ii < wpt.size() - 1; ++ii) {
	wpt[ii]->setNeighbors(wpt[ii-1], wpt[ii+1]);
      }
      wpt[wpt.size() - 1]->setNeighbors(wpt[wpt.size() - 2], goal);
    
      path_.push_back(start);
      for (size_t ii(0); ii < wpt.size(); ++ii) {
	path_.push_back(wpt[ii]);
      }
      path_.push_back(goal);
    
      for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
	(*ii)->init(state, Vector::Zero(state.size()));
      }
    }
      
      
    void InteractiveBlender::
    update()
    {
      z_angle_ = atan2(eestartori_.point_[1] - eestart_.point_[1], eestartori_.point_[0] - eestart_.point_[0]);
      Blender::update();
    }
    
    
    void InteractiveBlender::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      for (path_t::const_reverse_iterator ii(path_.rbegin()); ii != path_.rend(); ++ii) {
	CairoDrawable const * wpt(dynamic_cast<CairoDrawable const *>(*ii));
	if (wpt) {
	  wpt->draw(cr, weight, pixelsize);
	}
      }
	
      for (size_t ii(0); ii < handles_.size(); ++ii) {
	handles_[ii]->draw(cr, weight, pixelsize);
      }
	
      cairo_set_source_rgba(cr, 0.0, 1.0, 0.0, 0.5);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr, eestart_.point_[0], eestart_.point_[1]);
      cairo_line_to(cr, eestartori_.point_[0], eestartori_.point_[1]);
      cairo_stroke(cr);
    }
    
  }
  
}
