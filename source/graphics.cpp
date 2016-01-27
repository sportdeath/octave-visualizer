#include <cmath>
#include <SFML/Graphics.hpp>

#include "graphics.hpp"

void numToColor(float in, int * r, int * g, int * b) {
  *r = sin(2 * M_PI * in) * 127 + 128;
  *g = sin(2 * M_PI * in + 2*M_PI/3.) * 127 + 128;
  *b = sin(2 * M_PI * in + 4*M_PI/3.) * 127 + 128;
}

sf::Vector2f edgePosition(float cx, float cy, double ang) {

  if ( ang - 2 * M_PI >= atan2(-cy, cx) 
    || ang < atan2(cy, cx) ) {

    return sf::Vector2f( 2*cx, cy - cx * tan(ang) );

  } else if ( ang >= atan2(cy, cx)
      && ang < atan2(cy, -cx) ) {

    return sf::Vector2f(cx - cy * tan(ang - M_PI/2), 0);

  } else if ( ang >= atan2(cy, -cx)
      && ang - 2 * M_PI < atan2(-cy, -cx) ) {

    return sf::Vector2f(0, cy + cx * tan(ang - M_PI));

  } else if (ang - 2 * M_PI >= atan2(-cy, -cx)
      && ang - 2 * M_PI < atan2(-cy, cx) ) {

    return sf::Vector2f(cx + cy * tan(ang - 3*M_PI/2.) , 2*cy);

  }

  return sf::Vector2f(cx, cy);

}

sf::ConvexShape createCenteredTriangle(sf::Vector2u screen, double startAng, double widthAng) {
  sf::Vector2f center, right, left;

  center.x = screen.x/2.f;
  center.y = screen.y/2.f;

  left = edgePosition(center.x, center.y, fmod(startAng + widthAng, 2*M_PI));
  right = edgePosition(center.x, center.y, startAng);

  sf::ConvexShape tri;
  if (left.x != right.x && left.y != right.y) {

    tri = sf::ConvexShape(4);
    sf::Vector2f corner;

    if (right.x == screen.x) {
      corner.x = screen.x; corner.y = 0;
    } else if (right.y == 0) {
      corner.x = 0; corner.y = 0;
    } else if (right.x == 0) {
      corner.x = 0; corner.y = screen.y;
    } else if (right.y == screen.y) {
      corner.x = screen.x; corner.y = screen.y;
    }
  
    tri.setPoint(3, corner);
  } else {
    tri = sf::ConvexShape(3);
  }

  tri.setPoint(1,center);
  tri.setPoint(2,left);
  tri.setPoint(0,right);

  return tri;
}
