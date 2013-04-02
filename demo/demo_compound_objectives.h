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

#ifndef KINEMATIC_OBJECTIVES_DEMO_COMPOUND_OBJECTIVES_HPP
#define KINEMATIC_OBJECTIVES_DEMO_COMPOUND_OBJECTIVES_HPP

#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/joint_limit_objective.h>
#include <kinematic_objectives/obstacle_objective.h>
#include <kinematic_objectives/link_position_objective.h>
#include <kinematic_objectives/point_attraction_objective.h>
#include <kinematic_objectives/point_repulsion_objective.h>
#include <kinematic_objectives/joint_damping_objective.h>
#include "planar_orientation_objective.h"
#include "planar_robot.h"
#include "planar_distance.h"
#include "cairo_drawable.h"


namespace kinematic_objectives {
  
  namespace demo {
    
    class InteractiveBlender;
    
    
    class InteractionHandle
      : public CairoDrawable
    {
    public:
      InteractionHandle(double radius, double red, double green, double blue, double alpha);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      Vector point_;
      double radius_, red_, green_, blue_, alpha_;
    };
    
    
    class BaseCompoundObjective
      : public CompoundObjective,
	public CairoDrawable
    {
    public:
      BaseCompoundObjective(InteractiveBlender const & blender,
		   InteractionHandle const & repulsor,
		   double const & z_angle);
    
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct
      
      //// XXXX protected or so...
      
      PlanarRobot robot_; // XXXX keep this before any constraints so we can use its values for initializing them
      PlanarDistance distance_api_;
      
      InteractionHandle const & repulsor_;
      double const & z_angle_;
      
      JointLimitObjective joint_limits_;
    
      ObstacleObjective avoid_base_;
      ObstacleObjective avoid_ellbow_;
      ObstacleObjective avoid_wrist_;
      ObstacleObjective avoid_ee_;
    
      PlanarOrientationObjective orient_ee_;
    
      PointRepulsionObjective repulse_base_;
      PointRepulsionObjective repulse_ellbow_;
      PointRepulsionObjective repulse_wrist_;
      PointRepulsionObjective repulse_ee_;
    
      JointDampingObjective joint_damping_;
    };
    
    
    class NormalCompoundObjective
      : public BaseCompoundObjective
    {
    public:
      NormalCompoundObjective(InteractiveBlender const & blender,
		     InteractionHandle const & repulsor,
		     double const & z_angle);
      
      virtual ~NormalCompoundObjective();
      
      /**
	 \note Do not use in production code: calls exit() on error.
      */
      virtual void init(Vector const & position, Vector const & velocity);
      
      // virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct
      
      /**
	 \note Do not use in production code: calls exit() on error.
      */
      void setNeighbors(BaseCompoundObjective const * prev,
			BaseCompoundObjective const * next);
      
      //// XXXX protected or so...
      
      BaseCompoundObjective const * prev_;
      BaseCompoundObjective const * next_;
      
      vector<PointAttractionObjective*> attract_prev_;
      vector<PointAttractionObjective*> attract_next_;
    };
    
    
    class BoundaryCompoundObjective
      : public BaseCompoundObjective
    {
    public:
      BoundaryCompoundObjective(InteractiveBlender const & blender,
		       InteractionHandle const & repulsor,
		       double const & z_angle,
		       Vector const * eegoal,
		       Vector const * baseattractor);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void preUpdateHook(); // rfct

      //// XXXX protected or so...
      
      LinkPositionObjective eeobjective_;
      PointAttractionObjective attract_base_;
      Vector const * eegoal_;
      Vector const * baseattractor_;
    };
    
  }
  
}

#endif // KINEMATIC_OBJECTIVES_DEMO_COMPOUND_OBJECTIVES_HPP
