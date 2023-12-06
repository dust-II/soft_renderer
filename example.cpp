#include "rasterization.h"
#include <iostream>
int main() {
  raster::Scene scene;
  raster::Vec3f pos_raw{1, 1, 50};
  std::optional<raster::Vec2i> pos = scene.rasterize(pos_raw);
  std::cout << "******************************\n";
  std::cout << "3D space:\n";
  std::cout << "x:" << pos_raw.x << " y:" << pos_raw.y << " z:" << pos_raw.z
            << "\n";

  std::cout << "screen view space:\n";
  if (pos) {
    std::cout << "x=" << pos.value().x << "\n";
    std::cout << "y=" << pos.value().y << "\n";
  } else {
    std::cout << "point not in view of crema\n";
  }
}
