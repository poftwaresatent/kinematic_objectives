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

#include <kinematic_objectives/integrator.h>

namespace kinematic_objectives {
  
  Integrator::
  Integrator(double stepsize, Vector const & qd_max)
    : stepsize_(stepsize),
      qd_max_(qd_max)
  {
  }
  
  
  void Integrator::
  compute(Vector const & bias,
	  Vector const & q_in,
	  Vector const & qd_in,
	  Vector & q_out,
	  Vector & qd_out) const
  {
    static double const kp(100.0);
    static double const kd(20.0);
    
    Vector qdd(kp * bias);
    
    double saturation(0.0);
    for (Vector::Index ii(0); ii < qdd.size(); ++ii) {
      double const si(fabs((qdd[ii] / qd_max_[ii]) / kd));
      if (si > saturation) {
	saturation = si;
      }
    }
    if (saturation > 1.0) {
      qdd /= saturation;
    }
    
    qdd -= kd * qd_in;
    
    qd_out = qd_in + stepsize_ * qdd;
    q_out = q_in + stepsize_ * qd_out;
  }
  
}
