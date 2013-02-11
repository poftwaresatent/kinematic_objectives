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

#include "example_robot.hpp"
#include <gdk/gdk.h>
#include <err.h>


namespace kinematic_elastic {
  
  
  ExampleRobot::
  ExampleRobot()
    : radius_(0.5),
      len_a_(0.8),
      len_b_(0.6),
      len_c_(0.3),
      pos_a_(3),
      pos_b_(3),
      pos_c_(3)
  {
  }
  
  
  Vector const & ExampleRobot::
  getPosition() const
  {
    return position_;
  }
  
  
  Vector const & ExampleRobot::
  getVelocity() const
  {
    return velocity_;
  }
  
  
  Transform ExampleRobot::
  frame(size_t node) const
  {
    Transform tf(Transform::Identity());
    switch (node) {
    case 0:
      tf.translation() << position_[0], position_[1], 0.0;
      break;
    case 1:
      tf.translation() << position_[0], position_[1], 0.0;
      tf.linear() << c2_, -s2_, 0.0, s2_, c2_, 0.0, 0.0, 0.0, 1.0;
      break;
    case 2:
      tf.translation() << pos_a_[0], pos_a_[1], 0.0;
      tf.linear() << c23_, -s23_, 0.0, s23_, c23_, 0.0, 0.0, 0.0, 1.0;
      break;
    case 3:
      tf.translation() << pos_b_[0], pos_b_[1], 0.0;
      tf.linear() << c234_, -s234_, 0.0, s234_, c234_, 0.0, 0.0, 0.0, 1.0;
      break;
    default:
      errx (EXIT_FAILURE, "ExampleRobot::frame() called on invalid node %zu", node);
    }
    return tf;
  }
  
  
  Matrix ExampleRobot::
  computeJxo(size_t node, Vector const & gpoint) const
  {
    Matrix Jxo(Matrix::Zero(6, 5));
    switch (node) {
    case 3:
      Jxo(0, 4) = pos_b_[1] - gpoint[1];
      Jxo(1, 4) = gpoint[0] - pos_b_[0];
      Jxo(5, 4) = 1.0;
    case 2:
      Jxo(0, 3) = pos_a_[1] - gpoint[1];
      Jxo(1, 3) = gpoint[0] - pos_a_[0];
      Jxo(5, 3) = 1.0;
    case 1:
      Jxo(0, 2) = position_[1] - gpoint[1];
      Jxo(1, 2) = gpoint[0]    - position_[0];
      Jxo(5, 2) = 1.0;
    case 0:
      Jxo(0, 0) = 1.0;
      Jxo(1, 1) = 1.0;
      break;
    default:
      errx (EXIT_FAILURE, "Robot::computeJxo() called on invalid node %zu", node);
    }
    return Jxo;
  }
  
  
  void ExampleRobot::
  update(Vector const & position, Vector const & velocity)
  {
    if (position.size() != 5) {
      errx (EXIT_FAILURE, "Robot::update(): position has %zu DOF (but needs 5)", (size_t) position.size());
    }
    if (velocity.size() != 5) {
      errx (EXIT_FAILURE, "Robot::update(): velocity has %zu DOF (but needs 5)", (size_t) velocity.size());
    }
    position_ = position;
    velocity_ = velocity;
    
    c2_ = cos(position_[2]);
    s2_ = sin(position_[2]);
    ac2_ = len_a_ * c2_;
    as2_ = len_a_ * s2_;
    
    q23_ = position_[2] + position_[3];
    c23_ = cos(q23_);
    s23_ = sin(q23_);
    bc23_ = len_b_ * c23_;
    bs23_ = len_b_ * s23_;
    
    q234_ = q23_ + position_[4];
    c234_ = cos(q234_);
    s234_ = sin(q234_);
    cc234_ = len_c_ * c234_;
    cs234_ = len_c_ * s234_;
    
    pos_a_ <<
      position_[0] + ac2_,
      position_[1] + as2_,
      0.0;
    pos_b_ <<
      pos_a_[0] + bc23_,
      pos_a_[1] + bs23_,
      0.0;
    pos_c_ <<
      pos_b_[0] + cc234_,
      pos_b_[1] + cs234_,
      0.0;
  }
  
}
