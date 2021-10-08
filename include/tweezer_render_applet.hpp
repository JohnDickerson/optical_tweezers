#ifndef TWEEZER_RENDER_APPLET_HPP
#define TWEEZER_RENDER_APPLET_HPP

#include "defines.h"

#include <vlut/Applet.hpp>
#include <vlut/GeometryPrimitives.hpp>
#include <vl/SceneManagerActorTree.hpp>
#include <vl/Actor.hpp>
#include <vl/Effect.hpp>
#include <vl/Time.hpp>
#include <vl/Font.hpp>
#include <vl/Text.hpp>

#include "position_log.hpp"


class tweezer_render_applet : public vlut::Applet {
public:
  tweezer_render_applet(const vl::fvec2& grid_dim_x,
                        const vl::fvec2& grid_dim_y,
                        const vl::fvec2& grid_dim_z,
                        position_log& plog);
  void shutdown();
  void initEvent();
  void keyPressEvent(unsigned short, vl::EKey key);
  void runEvent();
  void run();

private:
  void generate_particles( vl::ref<vl::Light> light );
  void generate_axes( vl::ref<vl::Light> light );
  void generate_laser(vl::ref<vl::Light> light);

  void update_particle(const size_t& pid, const size_t& ts);
  void update_particles(const size_t& ts);


  vl::fvec2 _xbounds;
  vl::fvec2 _ybounds;
  vl::fvec2 _zbounds;

  position_log& _plog;

  vl::Collection< vl::Actor > _particles;
  vl::Time _timer;

  vl::ref<vl::Text> _title_text;
  vl::ref<vl::Text> _info_text;

  std::pair<size_t,size_t> _ts_counter;

};

#endif //TWEEZER_RENDER_APPLET_HPP
