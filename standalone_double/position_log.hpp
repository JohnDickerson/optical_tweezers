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

#include "constants.h"
#include "double3.h"


class position_log {
  friend std::ostream& operator<<(std::ostream& os, const position_log& plog);

public:
  position_log(const size_t& num_particles);
  position_log(const size_t& num_particles, const size_t& num_timesteps);
  ~position_log();

  double3S g;
  void append(const size_t& particle_id, double3S &pos);

	// Return
	double3S get_position(const size_t& particle_id);
  double3S get_position(const size_t& particle_id, const size_t& timestamp);

  size_t num_particles() const;
  size_t num_timesteps() const;

private:
  std::vector< std::vector<double3S> > _pos_log;
};


#endif // POSITION_LOG_HPP
