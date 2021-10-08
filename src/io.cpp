#include "io.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>

#include "constants.h"
//#include <cutil_math.h>

namespace tweezers
{
  namespace io
  {

    class forcemap_comp {
    public:
      bool operator() (const pair<float2,float2>& p1,
                       const pair<float2,float2>& p2) {
        if ( abs(p1.first.y - p2.first.y) < std::numeric_limits<float>::epsilon() ) {
          return p1.first.x < p2.first.x;
        } else {
          return p1.first.y < p2.first.y;
        }
      }

    };

    forcemap* load_force_data(std::string path, double dist_scale, double force_scale) {
		cerr << "# loading force data from : " << path << "\n";
      std::ifstream f_in;
      f_in.open(path.c_str(), std::ios::in);
			
      if(f_in.is_open())
        {
          forcemap* loaded_data = new forcemap();
          std::string line;
          double y,z,fy,fz;
          while( f_in >> y >> z >> fy >> fz ) {

	
            // Scale distances microns -> meters
            y *= dist_scale; z *= dist_scale;
            // Scale forces picoNewtons -> Newton
            // fy *= force_scale; fz *= force_scale;
            loaded_data->push_back( make_pair( make_float2(y,z), make_float2(fy,fz) ) );
            /*
            // Scale distances microns -> meters
            y *= dist_scale; z *= dist_scale;
					
            // We are keying on ints
            int key_y, key_z;
            key_y = round(y * constants::grid_inv_spacing);
            key_z = round(z * constants::grid_inv_spacing);
			
            // Scale forces picoNewtons -> Newton
            fy *= force_scale; fz *= force_scale;
					
            if( loaded_data->find(key_y) == loaded_data->end() ) {
            //	      if( !(*loaded_data)[key_y] )
            (*loaded_data)[key_y] = new map<int, forcepair>();
            }
			
            forcepair fp = {fy, fz};
            (*(*loaded_data)[key_y])[key_z] = fp;
            */
          }
			
          f_in.close();
          forcemap_comp fmcomp;
          std::sort( loaded_data->begin(), loaded_data->end(), fmcomp );
	  				
          return loaded_data;

        } else {
        std::cout << "Error opening file " << path << std::endl;
        return 0;
      }

    }
		
    void print_time_step(vector<particle*> *particles, laser* laserbeam, double time)
    {
      cout << endl << endl;
      cout << "# TIME (s)" << endl;
      cout << time << endl;
      cout << "# LASER" << endl;
      cout << *laserbeam << endl;
      cout << "# PARTICLES" << endl;
			
      for( vector<particle*>::iterator p_it = particles->begin(); p_it != particles->end(); p_it++)
	{
	  cout << **p_it << endl;
	}
    }
    /*		
    void print_loaded_force_data(forcemap* force_data)
    {
      forcemap::iterator y_it = force_data->begin();
      for(; y_it != force_data->end(); y_it++)
	{
	  map<int, forcepair>* y_row = y_it->second;
	  map<int, forcepair>::iterator z_it = y_row->begin();
	  for(; z_it != y_row->end(); z_it++)
	    {
	      forcepair fp = z_it->second;
	      cout << y_it->first << " " << z_it->first << " " << fp.fy << " " << fp.fz << endl;
	    }
	}
			
    }
    */	
    void delete_loaded_force_data(forcemap* force_data)
    {
      /*
	forcemap::iterator y_it = force_data->begin();
	for(; y_it != force_data->end(); y_it++)
	{
	delete y_it->second;
	}
      */
      delete force_data;
    }

    void load_config(string path, config_map& cm) {
      ifstream cm_in(path.c_str());
      if (cm_in.is_open()) {
        string key, value;
        while (cm_in >> key >> value) {
          cm[key] = value;
        }
      }
      cm_in.close();
    }
		

  }
}

