#include "io.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include "cutil_math.h"
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

void load_force_data(std::string path, double dist_scale, double force_scale, 
                     Array< blitz::TinyVector<double, 2>, 2 >& forcemap ) {

  forcemap.resize(113, 81);

  cerr << "loading force data from : " << path << "\n";
  std::ifstream f_in;
  f_in.open(path.c_str(), std::ios::in);

  if(f_in.is_open())
    {
      tforcemap loaded_data;
      //      forcemap* loaded_data = new forcemap();
      std::string line;
      double y,z,fy,fz;
      while( f_in >> y >> z >> fy >> fz ) {

	
        // Scale distances microns -> meters
        y *= dist_scale; z *= dist_scale;
        // Scale forces picoNewtons -> Newton
        // fy *= force_scale; fz *= force_scale;
        loaded_data.push_back( make_pair( make_float2(y,z), make_float2(fy,fz) ) );
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
      std::sort( loaded_data.begin(), loaded_data.end(), fmcomp );
      tforcemap::iterator it = loaded_data.begin();
	  	
      std::cerr << "num data points = " << loaded_data.size() << "\n";
      std::cerr << "expected data points = " << 81*113 << "\n";
      
      std::ofstream f_out;
      f_out.open("force_sorted_1.txt", std::ios::out);

      if(f_out.is_open()) {

          for (int y = 0; y < 113; ++y) {
            for (int z = 0; z < 81; ++z) {
              f_out << "F[" << y << ", " << z << "] = "
                    << static_cast<double>(it->second.x) << ", " 
                    << static_cast<double>(it->second.y) << " @ "
                    << static_cast<double>(it->first.x) << ", " 
                    << static_cast<double>(it->first.y) << "\n";
              ++it;
            }
          }
        
      }
      f_out.close();

      if (loaded_data.size() != 81*113) {
        std::exit(-1);
      }

      it = loaded_data.begin();
      for (int y = 0; y < 113; ++y) {
        for (int z = 0; z < 81; ++z) {
          forcemap(y,z)[0] = static_cast<double>(it->second.x);
          forcemap(y,z)[1] = static_cast<double>(it->second.y);
          ++it;
        }
      }

      // std::ofstream f_out2;
      // f_out2.open("force_sorted_2.txt", std::ios::out);

      // if(f_out2.is_open()) {
      //   float y_min = 0;
      //   float y_max = 20.0;
      //   float z_min = -20.0;
      //   float z_max = 8.0;

      //   for (int y = 0; y < 80; ++y) {
      //     for (int z = 0; z < 112; ++z) {
            
      //       forcemap(y,z)[0] = static_cast<double>(it->second.x);
      //       forcemap(y,z)[1] = static_cast<double>(it->second.y);
      //       ++it;
      //     }
      //   }

      // }
      // f_out2.close();

    } else {
    std::cout << "Error opening file " << path << std::endl;
  }

}
		
