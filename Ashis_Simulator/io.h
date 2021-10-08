#ifndef __TWEEZERS_IO_H__
#define __TWEEZERS_IO_H__

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <cutil_inline.h>


using namespace blitz;
using namespace std;

//typedef map<int, map<int, forcepair>*> forcemap;
typedef vector< pair<float2, float2> > tforcemap;
typedef map<string,string> config_map;

void load_force_data(string path, double dist_scale, double force_scale, 
                     Array< TinyVector<double, 2>, 2 >& forcemap );

#endif

