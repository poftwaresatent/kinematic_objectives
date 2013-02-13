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

#ifndef KINEMATIC_ELASTIC_EXAMPLE_WAYPOINTS_HPP
#define KINEMATIC_ELASTIC_EXAMPLE_WAYPOINTS_HPP

#include "waypoint.hpp"
#include "joint_limit_constraint.hpp"
#include "obstacle_constraint.hpp"
#include "position_control.hpp"
#include "point_attraction.hpp"
#include "point_repulsion.hpp"
#include "posture_damping.hpp"
#include "example_orientation_control.hpp"
#include "example_robot.hpp"
#include "example_distance_api.hpp"
#include "cairo_drawable.hpp"


namespace kinematic_elastic {
  
  namespace example {
    
    class InteractiveElastic;
    
    
    class InteractionHandle
      : public CairoDrawable
    {
    public:
      InteractionHandle(double radius, double red, double green, double blue, double alpha);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      Vector point_;
      double radius_, red_, green_, blue_, alpha_;
    };
    
    
    class BaseWaypoint
      : public Waypoint,
	public CairoDrawable
    {
    public:
      BaseWaypoint(InteractiveElastic const & elastic,
		   InteractionHandle const & repulsor,
		   double const & z_angle);
    
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct
      
      //// XXXX protected or so...
      
      PlanarRobot robot_; // XXXX keep this before any constraints so we can use its values for initializing them
      PlanarDistanceAPI distance_api_;
      
      InteractionHandle const & repulsor_;
      double const & z_angle_;
      
      JointLimitConstraint joint_limits_;
    
      ObstacleConstraint avoid_base_;
      ObstacleConstraint avoid_ellbow_;
      ObstacleConstraint avoid_wrist_;
      ObstacleConstraint avoid_ee_;
    
      OrientationControl orient_ee_;
    
      PointRepulsion repulse_base_;
      PointRepulsion repulse_ellbow_;
      PointRepulsion repulse_wrist_;
      PointRepulsion repulse_ee_;
    
      PostureDamping joint_damping_;
    };
    
    
    class NormalWaypoint
      : public BaseWaypoint
    {
    public:
      NormalWaypoint(InteractiveElastic const & elastic,
		     InteractionHandle const & repulsor,
		     double const & z_angle);
      
      virtual ~NormalWaypoint();
      
      /**
	 \note Do not use in production code: calls exit() on error.
      */
      virtual void init(Vector const & position, Vector const & velocity);
      
      // virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct
      
      /**
	 \note Do not use in production code: calls exit() on error.
      */
      void setNeighbors(BaseWaypoint const * prev,
			BaseWaypoint const * next);
      
      //// XXXX protected or so...
      
      BaseWaypoint const * prev_;
      BaseWaypoint const * next_;
      
      vector<PointAttraction*> attract_prev_;
      vector<PointAttraction*> attract_next_;
    };
    
    
    class BoundaryWaypoint
      : public BaseWaypoint
    {
    public:
      BoundaryWaypoint(InteractiveElastic const & elastic,
		       InteractionHandle const & repulsor,
		       double const & z_angle,
		       Vector const * eegoal,
		       Vector const * baseattractor);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct

      //// XXXX protected or so...
      
      PositionControl eetask_;
      PointAttraction attract_base_;
      Vector const * eegoal_;
      Vector const * baseattractor_;
    };
    
  }
  
}

#endif // KINEMATIC_ELASTIC_EXAMPLE_WAYPOINTS_HPP
