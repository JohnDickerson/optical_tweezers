#ifndef __util__hpp
#define __util__hpp

#include <boost/format.hpp>

namespace util{

  std::string to_s(const float3& v) {
    using boost::format;
    return str(format("<%1%, %2%, %3%>") % v.x % v.y % v.z);
  }

  std::string to_s(const float2& v) {
    using boost::format;
    return str(format("<%1%, %2%>") % v.x % v.y);
  }

  std::string to_s(const float4& v) {
    using boost::format;
    return str(format("<%1%, %2%, %3%, %4%>") % v.x % v.y % v.z % v.w);
  }

};

#endif
