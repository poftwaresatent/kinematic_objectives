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

#ifndef KINEMATIC_OBJECTIVES_DEMO_INTERACTIVE_COMPOUNDS_HPP
#define KINEMATIC_OBJECTIVES_DEMO_INTERACTIVE_COMPOUNDS_HPP

#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/joint_limit_objective.h>
#include <kinematic_objectives/obstacle_objective.h>
#include <kinematic_objectives/point_attraction_objective.h>
#include <kinematic_objectives/point_repulsion_objective.h>
#include "planar_orientation_objective.h"
#include "planar_robot.h"
#include "planar_distance.h"
#include "interaction_handle.h"


namespace kinematic_objectives {
  
  namespace demo {
    
    
    /**
       \todo [low] attributes should be protected or private
    */
    class InteractiveCompound
      : public CairoDrawable
    {
    public:
      explicit InteractiveCompound(PlanarRobot & robot);
    
      virtual void init(double gui_dimx, double gui_dimy);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void update();
      
      InteractionHandle h_ee_;
      InteractionHandle h_ee_ori_;
      InteractionHandle h_base_;
      InteractionHandle h_repulsor_;
      InteractionHandle h_obstacle_;
      
      vector<InteractionHandle*> handles_;
      
      PlanarRobot & robot_;
      PlanarDistance distance_api_;
      CompoundObjective compound_objective_;
      
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
    };
    
    
    /**
       \todo [low] attributes should be protected or private
    */
    class ElasticLinksCompound
      : public InteractiveCompound
    {
    public:
      explicit ElasticLinksCompound(PlanarRobot & robot);
      
      virtual void init(double gui_dimx, double gui_dimy);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void update();
      
      InteractionHandle h2_ee_;
      InteractionHandle h1_wrist_;
      InteractionHandle h2_wrist_;
      InteractionHandle h1_ellbow_;
      InteractionHandle h2_ellbow_;
      InteractionHandle h2_base_;
      
      PointAttractionObjective ee_left_;
      PointAttractionObjective ee_right_;
      PointAttractionObjective wrist_left_;
      PointAttractionObjective wrist_right_;
      PointAttractionObjective ellbow_left_;
      PointAttractionObjective ellbow_right_;
      PointAttractionObjective base_left_;
      PointAttractionObjective base_right_;
    };
    
    
    /**
       \todo [low] attributes should be protected or private
    */
    class EEGoalCompound
      : public InteractiveCompound
    {
    public:
      explicit EEGoalCompound(PlanarRobot & robot);
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      virtual void update();

      PointAttractionObjective attract_ee_;
      PointAttractionObjective attract_base_;
    };
    
  }
  
}

#endif // KINEMATIC_OBJECTIVES_DEMO_INTERACTIVE_COMPOUNDS_HPP
