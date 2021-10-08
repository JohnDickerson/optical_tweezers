#include <blitz/array.h>
#include <blitz/tinyvec.h>

namespace tweezers
{
  namespace interpolation
  {
    bool bi_interp_lookup ( const Array< blitz::TinyVector<double, 2>, 2 >& force_array,
                            double r, double z, TinyVector<double, 2>& interp_res ) {

      double rmin = 0.0e-6;
      double rmax = 20.0e-6;
      double zmin = -20.0e-6;
      double zmax = 8.0e-6;
    
      size_t nysamp = force_array.extent(0);
      size_t nzsamp = force_array.extent(1);


  
      double ind_r = ((r - rmin) / (rmax - rmin)) * (nysamp-1);
      double ind_z = ((z - zmin) / (zmax - zmin)) * (nzsamp-1);

      //int steps_per_micron = 4;    
      //double low_r = static_cast<double>(static_cast<int>(r * steps_per_micron)) / (steps_per_micron);
      //double low_z = static_cast<double>(static_cast<int>(z * steps_per_micron)) / (steps_per_micron);
  
      if (r > rmax || r < rmin || z > zmax || z < zmin) {
        interp_res[0] = 0.0;
        interp_res[1] = 0.0;
        return false;
      }

      int r_low = floor(ind_r);
      int r_hi = ceil(ind_r);
      int z_low = floor(ind_z);
      int z_hi = ceil(ind_z);

      double pr_low = ((static_cast<double>(r_low)/(nysamp-1)) * (rmax-rmin)) + rmin;
      double pr_hi = ((static_cast<double>(r_hi)/(nysamp-1)) * (rmax-rmin)) + rmin;
      double pz_low = ((static_cast<double>(z_low)/(nzsamp-1)) * (zmax-zmin)) + zmin;
      double pz_hi = ((static_cast<double>(z_hi)/(nzsamp-1)) * (zmax-zmin)) + zmin;

      // cerr << "\n\n" 
      //      << "pr_low = " << pr_low << " "
      //      << "pr_hi = " << pr_hi << " " 
      //      << "pz_low = " << pz_low << " " 
      //      << "pz_hi = " << pz_hi << "\n";

      //double ir_low = ((pr_low - rmin) / (rmax - rmin)) * (nysamp-1);
      //double ir_hi = ((pr_hi - rmin) / (rmax - rmin)) * (nysamp-1);
      //double iz_low = ((pz_low - zmin) / (zmax - zmin)) * (nzsamp-1);
      //double iz_hi = ((pz_hi - zmin) / (zmax - zmin)) * (nzsamp-1);

      // cerr << "ir_low = " << ir_low << " "
      //      << "ir_hi = " << ir_hi << " " 
      //      << "iz_low = " << iz_low << " " 
      //      << "iz_hi = " << iz_hi << "\n\n";


      TinyVector<double, 2> interp_rhi(0.0);
      TinyVector<double, 2> interp_rlow(0.0);
      TinyVector<double, 2> interp_z(0.0);
      TinyVector<double, 2> temp(0.0);

      double i0 = static_cast<double>(pr_hi - r) / (pr_hi - pr_low);
      double i1 = 1.0 - i0;
      if ( (r_hi == r_low) ) {
        i0 = 1.0;
        i1 = 0.0;
      }

      double j0 = static_cast<double>(pz_hi - z) / (pz_hi - pz_low);
      double j1 = 1.0 - j0;
      if ( (z_hi == z_low) ) {
        j0 = 1.0;
        j1 = 0.0;
      }

      // cerr << "ind_r = " << ind_r << ", ind_z = " << ind_z 
      //      << ", r_low = " << r_low << ", r_hi = " << r_hi 
      //      << ", z_low = " << z_low << ", z_hi = " << z_hi << "\n";
      // cerr << "i0 = " << i0 << ", i1 = " << i1 << ", j0 = " << j0 << ", j1 = " << j1 << "\n";


      interp_rlow = force_array(r_low, z_low);
      interp_rlow *= i0;
      temp = force_array(r_hi, z_low);
      temp *= i1;
      interp_rlow += temp;
      temp[0] = temp[1] = 0.0;

      // cerr << "interp_rlow = " << interp_rlow << " = "  
      //      << i0 << " *  f00 = " << force_array(z_low, r_low) << " "
      //      << i1  << " * f01 = " << force_array(z_low, r_hi) << "\n";
  
      interp_rhi = force_array(r_low, z_hi);
      interp_rhi *= i0;
      temp = force_array(r_hi, z_hi);
      temp *= i1;
      interp_rhi += temp;
      temp[0] = temp[1] = 0.0;

      // cerr << "interp_hi = " << interp_rhi << " = "  
      //      << i0 << " *  f10 = " << force_array(z_hi, r_low) << " "
      //      << i1  << " * f11 = " << force_array(z_hi, r_hi) << "\n";

      interp_z = interp_rlow;
      interp_z *= j0;
      temp = interp_rhi;
      temp *= j1;
      interp_z += temp;
      temp[0] = temp[1] = 0.0;

      interp_res[0] = interp_z[0] * 1e-12;
      interp_res[1] = interp_z[1] * 1e-12;

      //cerr << "[ " << r << ", " << z << "] ==> "
      //     << interp_res << "\n";


      return true;
    }

  }
}
