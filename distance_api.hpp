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

#ifndef KINEMATIC_ELASTIC_DISTANCE_API_HPP
#define KINEMATIC_ELASTIC_DISTANCE_API_HPP

#include "kinematic_elastic.hpp"


namespace kinematic_elastic {
  
  
  class DistanceAPI
  {
  public:
    struct reply_s {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      
      reply_s()
	: distance(-1.0), link_point(3), obstacle_point(3)
      {}
      
      /**
	 Distance between link and obstacle point, if available. Set
	 to -1.0 if there is no obstacle closer than
	 DistanceAPI::getMaxDistance().
      */
      double distance;
      
      /**
	 If a closest obstacle was found, then this is the point on
	 the link which lies (approximately) closest to the
	 obstacle. Note that DistanceAPI::handleRequest() yields
	 points with respect to the global frame.
      */
      Vector link_point;
      
      /**
	 If a closest obstacle was found, then this is the point on
	 the obstacle which lies (approximately) closest to the
	 link. Note that DistanceAPI::handleRequest() yields points
	 with respect to the global frame.
      */
      Vector obstacle_point;
    };
    
    virtual ~DistanceAPI() {}
    
    virtual double getMaxDistance() const = 0;
    virtual reply_s handleRequest(size_t link, double max_distance_hint) = 0;
  };
  
}

#endif // KINEMATIC_ELASTIC_DISTANCE_API_HPP
