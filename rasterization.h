#pragma once

#include "geometry.h"
#include <optional>

namespace raster {

inline constexpr bool isInside(const Vec4f &clip_pos) {
  float w = clip_pos.w;

  if (w == 0.0f)
    return false;
  if (clip_pos.z < 0.0f || clip_pos.z > w)
    return false;
  if (clip_pos.x < -w || clip_pos.x > w)
    return false;
  if (clip_pos.y < -w || clip_pos.y > w)
    return false;

  return true;
}

class Scene {
public:
  Scene() { mat_mvp_ = makeMVPMat(); }

  std::optional<Vec2i> rasterize(const Vec3f &point) {
    Vec4f pos{point.x, point.y, point.z, 1};
    Vec4f clip_pos = pos * mat_mvp_;

    if (!isInside(clip_pos)) {
      return {};
    }

    float rhw = 1.0f / clip_pos.w; // Reciprocal of the Homogeneous W
    Vec4f cvv_pos = clip_pos * rhw;

    Vec2f screen_pos_f{(cvv_pos.x + 1.0f) * viewport_width_ * 0.5f,
                       (1.0f - cvv_pos.y) * viewport_width_ * 0.5f};

    return Vec2i{static_cast<int32_t>(screen_pos_f.x + 0.5f),
                 static_cast<int32_t>(screen_pos_f.y + 0.5f)};
  }

private:
  Mat4x4f makeMVPMat() const {
    Mat4x4f mat_model =
        matRotate(0, 0, 0, 0) * matScale(1, 1, 1) * matTranslate(0, 0, 0);
    Mat4x4f mat_view = matLookat(eye_pos_, eye_at_, eye_up_);
    Mat4x4f mat_proj = matPerspecTive(
        fovy_, static_cast<float>(viewport_width_) / viewport_height_, near_,
        far_);
    return mat_model * mat_view * mat_proj;
  }

private:
  // camera
  const Vec3f eye_pos_{1, 1, 3}; // camera position
  const Vec3f eye_at_{0, 0, 0};  // camera direction
  const Vec3f eye_up_{0, 1, 0};  // camera up vector

  // project
  const float fovy_ = (3.1415926f / 180) * 45;
  const float near_ = 0.1f;
  const float far_ = 500.0f;

  // viewport
  const int32_t viewport_width_ = 800;
  const int32_t viewport_height_ = 480;

  Mat4x4f mat_mvp_;
};

} // namespace raster
