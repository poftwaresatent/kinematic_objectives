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

#include "planar_distance.h"
#include "planar_robot.h"
#include "interactive_blender.h"
#include <err.h>


namespace kinematic_objectives {
  
  namespace demo {
    
    PlanarDistance::
    PlanarDistance(PlanarRobot const & robot,
		      InteractiveBlender const & blender)
      : robot_(robot),
	blender_(blender)
    {
    }
     
    
    double PlanarDistance::
    computeMinimumSeparation(size_t link,
			     Vector & link_point,
			     Vector & obstacle_point)
      const
    {
      link_point.resize(3);
      obstacle_point = blender_.obstacle_.point_;
      Vector unit(3);
      double len;
      switch (link) {
      case 0:
	link_point << robot_.position_[0], robot_.position_[1], 0.0;
	unit = obstacle_point - link_point;
	len = unit.norm();
	unit /= len;
	link_point += robot_.radius_ * unit;
	obstacle_point -= blender_.obstacle_.radius_ * unit;
	return len - robot_.radius_ - blender_.obstacle_.radius_;
      case 1:
	link_point << robot_.position_[0], robot_.position_[1], 0.0;
	unit << robot_.c2_, robot_.s2_, 0.0;
	len = robot_.len_a_;
	break;
      case 2:
	link_point << robot_.pos_a_[0], robot_.pos_a_[1], 0.0;
	unit << robot_.c23_, robot_.s23_, 0.0;
	len = robot_.len_b_;
	break;
      case 3:
	link_point << robot_.pos_b_[0], robot_.pos_b_[1], 0.0;
	unit << robot_.c234_, robot_.s234_, 0.0;
	len = robot_.len_c_;
	break;
      default:
	errx (EXIT_FAILURE, "PlanarDistance::handleRequest() called on invalid link %zu", link);
      }

      Vector delta(obstacle_point - link_point);
      double const perp(delta.transpose() * unit);
      if (perp > 0.0) {
	if (perp >= len) {
	  link_point += unit * len;
	}
	else {
	  link_point += unit * perp;
	}
	delta = obstacle_point - link_point;
      }
      
      len = delta.norm();
      unit = delta / len;
      obstacle_point -= blender_.obstacle_.radius_ * unit;
      return len - blender_.obstacle_.radius_;
    }
    
  }

}
