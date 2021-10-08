#include "tweezer_render_applet.hpp"

#include <sstream>

#include <vl/VisualizationLibrary.hpp>
#include <vl/Effect.hpp>
#include <vl/FontManager.hpp>

tweezer_render_applet::tweezer_render_applet( const vl::fvec2& grid_dim_x,
                                              const vl::fvec2& grid_dim_y,
                                              const vl::fvec2& grid_dim_z,
                                              position_log& plog) :
  _xbounds(grid_dim_x), _ybounds(grid_dim_y), _zbounds(grid_dim_z),
  _plog(plog), _ts_counter( make_pair(0,plog.num_timesteps()-1) )
{}

void tweezer_render_applet::shutdown() {}
/*  
void tweezer_render_applet::set_mesh(vl::ref<vl::Geometry>& geom) {
  vl::ref<vl::Transform> xform = new vl::Transform;
  vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());
  
  vl::ref<vl::Light> light = new vl::Light(0);
 
  vl::ref<vl::Effect> fx = new vl::Effect;
  fx->shader()->enable(vl::EN_DEPTH_TEST);
  fx->shader()->gocMaterial()->setColorMaterialEnabled(true);
  fx->shader()->enable(vl::EN_LIGHTING);
  fx->shader()->setRenderState( light.get() );
  fx->setRenderRank(1);

  vl::ref<vl::Actor> act;
  act = sceneManager()->tree()->addActor( _mesh.get(), fx.get(), xform.get());

  trackball()->adjustView( vl::VisualizationLibrary::rendering()->as<vl::Rendering>(), 
			   vl::vec3(0,0,1), vl::vec3(0,1,0), 1.0f );  
  sceneManager()->computeBounds();
}
*/

void tweezer_render_applet::generate_laser(vl::ref<vl::Light> light) {
  vl::vec3 origin(0.0f, 0.0f, 0.0f);
  vl::vec3 xaxis(1.0, 0.0, 0.0);
  vl::vec3 yaxis(0.0, 1.0, 0.0); 
  vl::vec3 zaxis(0.0, 0.0, 1.0);

  vl::fvec4 color(0.85f, 0.1f, 0.1f, 0.7f);

  // Create the laser
  {
    vl::ref<vl::Geometry> laser = vlut::makeCylinder(origin, 0.25f, 23.0f);
    laser->computeNormals();

    vl::ref<vl::Geometry> laser_focus = vlut::makeCone(origin, 4.0f, 3.0f,
                                                       10, false);
    laser_focus->computeNormals();

    vl::mat4 tmat;
    tmat.setIdentity();
    tmat.rotate(90.0f, yaxis);
    tmat.translate(0.0f, 5.0f, 0.0f);

    vl::ref<vl::Transform> xform = new vl::Transform;
    xform->setLocalMatrix(tmat);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());

    vl::mat4 tmat_focus;
    tmat_focus.setIdentity();
    tmat_focus.rotate(180.0f, zaxis);
    tmat_focus.translate(0.0f, 0.0f, 0.0f);

    vl::ref<vl::Transform> xform_focus = new vl::Transform;
    xform_focus->setLocalMatrix(tmat_focus);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform_focus.get());

    vl::ref<vl::Effect> fx = new vl::Effect;
    fx->shader()->enable(vl::EN_DEPTH_TEST);
    fx->shader()->gocMaterial()->setDiffuse( color );
    fx->shader()->enable(vl::EN_LIGHTING);
    fx->shader()->setRenderState( light.get() );
    fx->shader()->enable(vl::EN_BLEND);

    vl::ref<vl::Actor> act;
    act = sceneManager()->tree()->addActor( laser.get(), fx.get(), xform.get() );
    act = sceneManager()->tree()->addActor( laser_focus.get(), fx.get(), xform_focus.get() );
  }

}

void tweezer_render_applet::generate_axes(vl::ref<vl::Light> light) {
  vl::vec3 origin(0.0f, 0.0f, 0.0f);
  vl::vec3 xaxis(1.0, 0.0, 0.0);
  vl::vec3 yaxis(0.0, 1.0, 0.0); 
  vl::vec3 zaxis(0.0, 0.0, 1.0);

  vl::fvec4 color(0.9f, 0.9f, 0.9f, 1.0f);

  vl::ref<vl::Texture> grid_cell = new vl::Texture("textures/grid_cell.png", vl::TF_RGBA, true, false);
  vl::TexParameter* params = grid_cell->getTexParameter();
  params->setGenerateMipmap(true);
  params->setMinFilter(vl::TPF_LINEAR_MIPMAP_LINEAR);
  params->setMagFilter(vl::TPF_LINEAR_MIPMAP_LINEAR);
  // Create the axes
  // XY Plane
  {
    vl::ref<vl::Geometry> xy_plane = vlut::makeGrid(origin, 
                                                    20.0f, 20.0f, 
                                                    100, 100,
                                                    true, 
                                                    vl::fvec2(0,0),vl::fvec2(20,20));
    xy_plane->computeNormals();

    vl::mat4 tmat;
    tmat.setIdentity();
    tmat.translate(0.0f, -5.0f, 0.0f);

    vl::ref<vl::Transform> xform = new vl::Transform;
    xform->setLocalMatrix(tmat);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());

    vl::ref<vl::Effect> fx = new vl::Effect;
    fx->shader()->enable(vl::EN_DEPTH_TEST);
    fx->shader()->gocMaterial()->setDiffuse( color );
    fx->shader()->enable(vl::EN_LIGHTING);
    fx->shader()->enable(vl::EN_CULL_FACE);
    fx->shader()->setRenderState( light.get() );
    fx->shader()->gocTextureUnit(0)->setTexture( grid_cell.get() );

    vl::ref<vl::Actor> act;
    act = sceneManager()->tree()->addActor( xy_plane.get(), fx.get(), xform.get() );
  }

  // YZ Plane
  {
    vl::ref<vl::Geometry> xy_plane = vlut::makeGrid(origin, 
                                                    23.0f, 20.0f, 
                                                    100, 100,
                                                    true,
                                                    vl::fvec2(0,0),vl::fvec2(23,20));

    xy_plane->computeNormals();

    vl::mat4 tmat;
    tmat.setIdentity();

    tmat.rotate(90.0f, zaxis);
    tmat.translate(10.0f, 6.5f, 0.0f);

    vl::ref<vl::Transform> xform = new vl::Transform;
    xform->setLocalMatrix(tmat);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());

    vl::ref<vl::Effect> fx = new vl::Effect;
    fx->shader()->enable(vl::EN_DEPTH_TEST);
    fx->shader()->gocMaterial()->setDiffuse( color );
    fx->shader()->enable(vl::EN_LIGHTING);
    fx->shader()->enable(vl::EN_CULL_FACE);
    fx->shader()->setRenderState( light.get() );
    fx->shader()->gocTextureUnit(0)->setTexture( grid_cell.get() );

    vl::ref<vl::Actor> act;
    act = sceneManager()->tree()->addActor( xy_plane.get(), fx.get(), xform.get() );
  }

  // XZ Plane
  {
    vl::ref<vl::Geometry> xy_plane = vlut::makeGrid(origin, 
                                                    20.0f, 23.0f, 
                                                    100, 100,
                                                    true,
                                                    vl::fvec2(0,0),vl::fvec2(20,23));

    xy_plane->computeNormals();

    vl::mat4 tmat;
    tmat.setIdentity();

    tmat.rotate(90.0f, xaxis);
    tmat.translate(0.0f, 6.5f, -10.0f);

    vl::ref<vl::Transform> xform = new vl::Transform;
    xform->setLocalMatrix(tmat);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());

    vl::ref<vl::Effect> fx = new vl::Effect;
    fx->shader()->enable(vl::EN_DEPTH_TEST);
    fx->shader()->gocMaterial()->setDiffuse( color );
    fx->shader()->enable(vl::EN_LIGHTING);
    fx->shader()->enable(vl::EN_CULL_FACE);
    fx->shader()->setRenderState( light.get() );
    fx->shader()->gocTextureUnit(0)->setTexture( grid_cell.get() );

    vl::ref<vl::Actor> act;
    act = sceneManager()->tree()->addActor( xy_plane.get(), fx.get(), xform.get() );
  }


}

void tweezer_render_applet::generate_particles(vl::ref<vl::Light> light) {

  vl::vec3 origin(0.0f, 0.0f, 0.0f);
  // Create the initial particles

  for (size_t pid = 0; pid < _plog.num_particles(); ++pid) {
   
    float3 pos = _plog.get_position(pid, _ts_counter.first);
    vl::mat4 tmat;
    tmat.setIdentity();
    tmat.translate(1e6*pos.x, -1e6*pos.z, 1e6*pos.y);

    vl::ref<vl::Transform> xform = new vl::Transform;
    xform->setLocalMatrix(tmat);
    vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->transform()->addChild(xform.get());

    vl::fvec4 color(rand()/static_cast<float>(RAND_MAX), rand()/static_cast<float>(RAND_MAX),
		   rand()/static_cast<float>(RAND_MAX), 1.0f);

    vl::ref<vl::Effect> fx = new vl::Effect;
    fx->shader()->enable(vl::EN_DEPTH_TEST);
    fx->shader()->gocMaterial()->setDiffuse( color );
    fx->shader()->enable(vl::EN_LIGHTING);
    fx->shader()->setRenderState( light.get() );

    vl::ref<vl::Geometry> psphere = 
      vlut::makeUVSphere(origin, 5.0, 20, 20);
    psphere->computeNormals();

    vl::ref<vl::Actor> act;
    act = sceneManager()->tree()->addActor( psphere.get(), fx.get(), xform.get() );
    _particles.push_back(act.get());
  }
  ++_ts_counter.first;
}

void tweezer_render_applet::initEvent() {
  _title_text = new vl::Text();
  _title_text->setDisplayListEnabled(true);
  _title_text->setFont( vl::VisualizationLibrary::fontManager()->acquireFont("font/bitstream-vera/Vera.ttf", 14) );
  _title_text->setMargin(5);
  _title_text->setViewportAlignment(vl::AlignTop | vl::AlignHCenter);
  _title_text->setAlignment(vl::AlignTop | vl::AlignHCenter);
  _title_text->setTextAlignment(vl::TextAlignCenter);
  _title_text->setColor(vlut::white);
  _title_text->setBackgroundColor(vl::fvec4(0,0,0,.75f));
  _title_text->setBackgroundEnabled(true);
  _title_text->setText("Visualization: 5 micron particles");
  _title_text->setDisplayListDirty(true);

  _info_text = new vl::Text();
  _info_text->setDisplayListEnabled(true);
  _info_text->setFont( vl::VisualizationLibrary::fontManager()->acquireFont("font/bitstream-vera/Vera.ttf", 14) );
  _info_text->setMargin(5);
  _info_text->setViewportAlignment(vl::AlignBottom | vl::AlignLeft);
  _info_text->setAlignment(vl::AlignBottom | vl::AlignLeft);
  _info_text->setTextAlignment(vl::TextAlignCenter);
  _info_text->setColor(vlut::white);
  _info_text->setBackgroundColor(vl::fvec4(0,0,0,.75f));
  _info_text->setBackgroundEnabled(true);
  std::stringstream ss(std::stringstream::in | std::stringstream::out);
  ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
  _info_text->setText(ss.str().c_str());
  _info_text->setDisplayListDirty(true);

  vl::ref< vl::Effect > text_fx = new vl::Effect;
  text_fx->shader()->enable(vl::EN_BLEND);
  text_fx->shader()->disable(vl::EN_DEPTH_TEST);
  text_fx->setRenderRank(10);

  vl::ref<vl::Actor> ttext_actor = new vl::Actor(_title_text.get(), text_fx.get()) ;
  sceneManager()->tree()->addActor(ttext_actor.get());
  vl::ref<vl::Actor> itext_actor = new vl::Actor(_info_text.get(), text_fx.get()) ;
  sceneManager()->tree()->addActor(itext_actor.get());

  vl::vec3 minvec(_xbounds[0],_zbounds[0],_ybounds[0]);
  vl::vec3 maxvec(_xbounds[1],_zbounds[1],_ybounds[1]);
  vl::AABB bb(minvec,maxvec);

  vl::vec3 middle = (minvec + maxvec) / 2.0;
  trackball()->setPivot(middle);
  trackball()->setTransform(NULL);


  //  trackball()->adjustView( bb, vl::vec3(0,0,1), vl::vec3(0,1,0), 1.0f );  
  sceneManager()->setBoundingBox(bb);

  vl::ref<vl::Light> light = new vl::Light(0);
  generate_laser(light);
  generate_axes(light);
  generate_particles(light);
}

void tweezer_render_applet::keyPressEvent(unsigned short, vl::EKey key) {
  switch(key) {

  case vl::Key_Up :
    if (_ts_counter.first < _ts_counter.second) {
      ++_ts_counter.first;
      update_particles(_ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }

    break;
  case vl::Key_Down :
    if (_ts_counter.first > 0) {
      --_ts_counter.first;
      update_particles(_ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }

    break;

  case vl::Key_P :
    if (_timer.isStarted(0)) {
      _timer.stop(0);
    } else {
      _timer.start(0);
    }
    break;

  case vl::Key_PageDown :
    {
      _ts_counter.first = 0;
      update_particles(_ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }
    break;


  case vl::Key_PageUp :
    {
      _ts_counter.first = _ts_counter.second;
      update_particles(_ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }
    break;


    /*
  case vl::Key_PageDown :
    if (_ts_counter.first > 0) {
      _ts_counter.first = std::max( static_cast<int>(0),
                                    static_cast<int>(_ts_counter.first) - 20 );
      igss::vis::update_colors(_mesh, *_igss_mat, _ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "IGSS timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }
    break;

  case vl::Key_PageUp :
    if (_ts_counter.first < _ts_counter.second-1) {
      _ts_counter.first = std::min( static_cast<int>(_ts_counter.second-1),
                                    static_cast<int>(_ts_counter.first + 20) );
      igss::vis::update_colors(_mesh, *_igss_mat, _ts_counter.first);
      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "IGSS timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);
    }
    break;
    */
  }
}

void tweezer_render_applet::runEvent() {
  vlut::Applet::runEvent();
}

void tweezer_render_applet::update_particle(const size_t& pid, const size_t& ts) {
  float3 pos = _plog.get_position(pid, _ts_counter.first);
  vl::mat4 tm; tm.setIdentity();
  tm.translate(1e6*pos.x, -1e6*pos.z, 1e6*pos.y);
  _particles[pid]->transform()->setLocalMatrix(tm);
}

void tweezer_render_applet::update_particles(const size_t& ts) {
  for (size_t pid = 0; pid < _particles.size(); ++pid) {
    update_particle(pid, ts);
  }
}

void tweezer_render_applet::run() {
  // The default action is good for now
  if (_timer.isStarted() && _timer.elapsed() > 0.032 ) {
    _timer.stop(0); // Later we'll do this in a less hacky way

    if (_ts_counter.first < _ts_counter.second) {
      ++_ts_counter.first;
      update_particles(_ts_counter.first);

      std::stringstream ss(std::stringstream::in | std::stringstream::out);
      ss << "Timestep " << _ts_counter.first << " / " << _ts_counter.second;
      _info_text->setText(ss.str().c_str());
      _info_text->setDisplayListDirty(true);

      _timer.start(0);// Later we'll do this in a less hacky way
    } 
    // Don't restart the timer if we've reached the last frame
  }

}


