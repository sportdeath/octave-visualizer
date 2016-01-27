#ifndef GRAPHICS
#define GRAPHICS

void numToColor(float in, int * r, int * g, int * b);

sf::Vector2f edgeCoordinates(float cx, float cy, double ang);

sf::ConvexShape createCenteredTriangle(sf::Vector2u screen, double startAng, double widthAng);

#endif
