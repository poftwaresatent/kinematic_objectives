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

#ifndef KINEMATIC_OBJECTIVES_DISTANCE_MODEL_HPP
#define KINEMATIC_OBJECTIVES_DISTANCE_MODEL_HPP

#include <kinematic_objectives/kinematic_objectives.h>


namespace kinematic_objectives {
  
  
  /**
     Functionality adapter that provides independence wrt the way
     distances are computed.
   */
  class DistanceModel
  {
  public:
    virtual ~DistanceModel() {}
    
    /**
       \return The minimum separation distance of the given link. This
       is a positive number (corresponding to the minimum distance) in
       case the link is collision-free, and a negative number in case
       of collision (with magnitude equal to the smallest displacement
       necessary to separate the link from the collision).
       
       \todo Does the moveit collision/distance API provide minimum
       separation distances?
       
       \param[out] link_point point on the link, in global
       coordinates, that is closest to the obstacle (in case of
       positive return value) or which would have to move the minimum
       amount (in case of negative return value). In case of multiple
       such points, a (hopefully) reasonable one will be selected.
       
       \param[out] obstacle_point point on the obstacle, in global
       coordinates, that is the returned distance away from \p
       link_point. As for the latter, a reasonable choice is made when
       there is more than one such point.
    */
    virtual double computeMinimumSeparation(size_t link,
					    Vector & link_point,
					    Vector & obstacle_point) const = 0;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_DISTANCE_MODEL_HPP
