#include "position_log.hpp"

position_log::position_log(const size_t& num_particles) :
  _pos_log( std::vector< std::vector<float3> >(num_particles) ) {}

position_log::position_log(const size_t& num_particles, const size_t& num_timesteps) :
  _pos_log( std::vector< std::vector<float3> >(num_particles) ) 
{
  typedef std::vector< std::vector<float3> >::iterator plog_it_t;
  for (plog_it_t it = _pos_log.begin(); it != _pos_log.end(); ++it) {
    it->reserve(num_timesteps);
  }
}

position_log::~position_log() {}

void position_log::append(const size_t& particle_id, float3& pos) {
  if (particle_id < _pos_log.size()) {
    _pos_log[particle_id].push_back(pos);
  }
}

float3 position_log::get_position(const size_t& particle_id)
{
	return get_position(particle_id, num_timesteps() - 1);
}


float3 position_log::get_position(const size_t& particle_id, const size_t& timestamp) {
  return _pos_log[particle_id][timestamp];
}

size_t position_log::num_particles() const {
  return _pos_log.size();
}

size_t position_log::num_timesteps() const {
  return _pos_log.front().size();
}

std::ostream& operator<<(std::ostream& os, const position_log& plog) {
  os << "@position log :: #particles : " << plog.num_particles() << " "
     << "#timesteps : " << plog.num_timesteps() << "\n";

  for (size_t ts = 0; ts < plog.num_timesteps(); ++ts) {
    os << "@timestep : " << ts << "\n";
    for (size_t pid = 0; pid < plog.num_particles(); ++pid ) {
      os << "@particle " << pid << " : (" 
         << plog._pos_log[pid][ts].x << ", "
         << plog._pos_log[pid][ts].y << ", "
         << plog._pos_log[pid][ts].z << ")\n";
    }
  }

  return os;
}
