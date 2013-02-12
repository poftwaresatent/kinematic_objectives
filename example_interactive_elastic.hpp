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

#ifndef KINEMATIC_ELASTIC_EXAMPLE_INTERACTIVE_ELASTIC_HPP
#define KINEMATIC_ELASTIC_EXAMPLE_INTERACTIVE_ELASTIC_HPP

#include "elastic.hpp"
#include "example_waypoints.hpp"


namespace kinematic_elastic {

  namespace example {
    
    class InteractiveElastic
      : public Elastic,
	public CairoDrawable
    {
    public:
      InteractiveElastic(double timestep,
			 ostream * dbgos,
			 string const & dbgpre);
      
      void init(Vector const & state);
      
      virtual void update();
      
      virtual void draw(cairo_t * cr, double weight, double pixelsize) const;
      
      //// XXXX protected or whatnot...
      
      InteractionHandle eestart_;
      InteractionHandle eestartori_;
      double z_angle_;
      InteractionHandle basestart_;
      InteractionHandle eegoal_;
      InteractionHandle basegoal_;
      InteractionHandle repulsor_;
      InteractionHandle obstacle_;
      
      vector<InteractionHandle*> handles_;
    };
    
  }

}

#endif // KINEMATIC_ELASTIC_EXAMPLE_INTERACTIVE_ELASTIC_HPP
