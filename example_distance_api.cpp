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

#include "example_distance_api.hpp"
#include "example_robot.hpp"
#include "example_interactive_elastic.hpp"
#include <err.h>


namespace kinematic_elastic {
  
  namespace example {
    
    PlanarDistanceAPI::
    PlanarDistanceAPI(PlanarRobot const & robot,
		      InteractiveElastic const & elastic)
      : robot_(robot),
	elastic_(elastic)
    {
    }
     
    
    double PlanarDistanceAPI::
    getMaxDistance() const
    {
      return numeric_limits<double>::max();
    }
    
    
    DistanceAPI::reply_s PlanarDistanceAPI::
    handleRequest(size_t link, double max_distance_hint)
    {
      DistanceAPI::reply_s reply;
      reply.obstacle_point = elastic_.obstacle_.point_;
      Vector unit(3);
      double len;
      switch (link) {
      case 0:
	reply.link_point << robot_.position_[0], robot_.position_[1], 0.0;
	reply.distance = (reply.obstacle_point - reply.link_point).norm();
	return reply;
      case 1:
	reply.link_point << robot_.position_[0], robot_.position_[1], 0.0;
	unit << robot_.c2_, robot_.s2_, 0.0;
	len = robot_.len_a_;
	break;
      case 2:
	reply.link_point << robot_.pos_a_[0], robot_.pos_a_[1], 0.0;
	unit << robot_.c23_, robot_.s23_, 0.0;
	len = robot_.len_b_;
	break;
      case 3:
	reply.link_point << robot_.pos_b_[0], robot_.pos_b_[1], 0.0;
	unit << robot_.c234_, robot_.s234_, 0.0;
	len = robot_.len_c_;
	break;
      default:
	errx (EXIT_FAILURE, "PlanarDistanceAPI::handleRequest() called on invalid link %zu", link);
      }

      Vector const delta(reply.obstacle_point - reply.link_point);
      double const perp(delta.transpose() * unit);
      if (perp <= 0.0) {
	reply.distance = delta.norm();
	return reply;
      }
      if (perp >= len) {
	reply.link_point += unit * len;
      }
      else {
	reply.link_point += unit * perp;
      }
      reply.distance = (reply.obstacle_point - reply.link_point).norm();
      return reply;
    }
    
  }

}
