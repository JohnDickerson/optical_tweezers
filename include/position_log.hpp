#ifndef POSITION_LOG_HPP
#define POSITION_LOG_HPP

/** 
 * Stores the positions of particles whose paths have
 * been calculated by our integrator.  Right now this 
 * class is just "stupid" and holds all the results in 
 * memory.  In the future it should do fancy things like
 * store results to disk and only buffer the most relevant 
 * positions, or even store the CUDA pointers to the most
 * current particle positions.
 */
#include <iostream>
#include <vector>

#include <cutil_inline.h>
#include "cutil_math.h"
#include "vector3d.h"

class position_log {
  friend std::ostream& operator<<(std::ostream& os, const position_log& plog);

public:
  position_log(const size_t& num_particles);
  position_log(const size_t& num_particles, const size_t& num_timesteps);
  ~position_log();

  void append(const size_t& particle_id, float3& pos);

	// Return
	float3 get_position(const size_t& particle_id);
  float3 get_position(const size_t& particle_id, const size_t& timestamp);

  size_t num_particles() const;
  size_t num_timesteps() const;

private:
  std::vector< std::vector<float3> > _pos_log;
};


#endif // POSITION_LOG_HPP
