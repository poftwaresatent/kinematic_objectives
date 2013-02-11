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

#ifndef KINEMATIC_ELASTIC_EXAMPLE_ROBOT_HPP
#define KINEMATIC_ELASTIC_EXAMPLE_ROBOT_HPP

#include "model.hpp"
#include <cairo/cairo.h>


namespace kinematic_elastic {
  
  
  class ExampleRobot
    : public Model
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
    ExampleRobot();
    
    virtual Vector const & getPosition() const;
    virtual Vector const & getVelocity() const;
    
    /**
       \note Do not use in production code: this method calls exit()
       when you specify an invalid node!
    */
    virtual Transform frame(size_t node) const;
  
    /**
       \note Do not use in production code: this method calls exit()
       when you specify an invalid node!
    */
    virtual Matrix computeJxo(size_t node, Vector const & gpoint) const;
    
    /**
       \note Do not use in production code: this method calls exit()
       when you specify an invalid position or velocity!
    */
    virtual void update(Vector const & position, Vector const & velocity);
    
    void draw(cairo_t * cr, double weight, double pixelsize) const;
    
    
    //// XXXX make these protected or whatnot...
  
    double const radius_;
    double const len_a_;
    double const len_b_;
    double const len_c_;
  
    Vector position_;
    Vector velocity_;
    Vector pos_a_;
    Vector pos_b_;
    Vector pos_c_;
  
    double c2_;
    double s2_;
    double c23_;
    double s23_;
    double c234_;
    double s234_;
    double q23_;
    double q234_;
    double ac2_;
    double as2_;
    double bc23_;
    double bs23_;
    double cc234_;
    double cs234_;
  };
  
}

#endif // KINEMATIC_ELASTIC_EXAMPLE_ROBOT_HPP
