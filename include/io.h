#ifndef __TWEEZERS_IO_H__
#define __TWEEZERS_IO_H__

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <blitz/array.h>
#include <blitz/tinyvec.h>

#include "laser.h"
#include "particle.h"

using namespace blitz;
using namespace std;

namespace tweezers
{
  namespace io
  {
    // Empirical force data for y, z coordinate
    typedef struct
    {
      double fy;
      double fz;
    } forcepair;
		
    //typedef map<int, map<int, forcepair>*> forcemap;
    typedef vector< pair<float2, float2> > forcemap;
    typedef map<string,string> config_map;

    forcemap* load_force_data(string path, double dist_scale, double force_scale);
    void load_config(string path, config_map& cm);

    void print_loaded_force_data(forcemap* force_data);
    void print_time_step(vector<particle*> *particles, laser* laserbeam, double time);
		
    void delete_loaded_force_data(forcemap* force_data);
  }
}

#endif

