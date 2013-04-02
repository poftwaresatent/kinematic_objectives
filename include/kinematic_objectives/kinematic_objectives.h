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

#ifndef KINEMATIC_OBJECTIVES_HPP
#define KINEMATIC_OBJECTIVES_HPP

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Geometry>
#include <cmath>


namespace kinematic_objectives {
  
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::Isometry3d Transform;
  
  using namespace std;
  
  void stackVector (Vector const & v1,
		    Vector const & v2,
		    Vector & vv);
  
  void stackMatrix (Matrix const & m1,
		    Matrix const & m2,
		    Matrix & mm);
  
  template<typename value_t>
  inline value_t bound(value_t lower, value_t value, value_t upper)
  {
    if (value < lower) {
      value = lower;
    }
    else if (value > upper) {
      value = upper;
    }
    return value;
  }
  
  
  static inline double normangle(double phi)
  {
    phi = fmod(phi, 2.0 * M_PI);
    if (phi > M_PI) {
      phi -= 2 * M_PI;
    }
    else if (phi < -M_PI) {
      phi += 2 * M_PI;
    }
    return phi;
  }
  
  static double const deg(M_PI / 180.);
  
}

#endif // KINEMATIC_OBJECTIVES_HPP
