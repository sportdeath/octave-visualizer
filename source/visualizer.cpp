#include <iostream>
#include <cmath>
#include <array>

#include <SFML/Graphics.hpp>

#include "visualizer.hpp"
#include "graphics.hpp"

Visualizer::Visualizer() {
  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  sf::ContextSettings settings;
  settings.antialiasingLevel = 4;

  window.create(
      sf::VideoMode(600, 600),
      "Visualizer",
      sf::Style::Titlebar |
      sf::Style::Resize |
      sf::Style::Close,
      settings);

  window.setVerticalSyncEnabled(true);

  norm = 0;
  prevMax = 0;
}

void Visualizer::setSpec(int index, double val) {
  spectrogram[index] = val;
}

double Visualizer::getSpec(int index) {
  return spectrogram[index];
}

void Visualizer::updateWindow() {
  //window.clear(sf::Color::Black);

  double delta = 2 * M_PI/SPEC_SIZE;
  int r, g, b;

  for (int i = 0; i < SPEC_SIZE; i++) {
    sf::ConvexShape tri = createCenteredTriangle(window.getSize(), delta*i, delta);

    double a = spectrogram[i];
    numToColor(i/(double)SPEC_SIZE, &r, &g, &b);
    tri.setFillColor(sf::Color(r*a,g*a,b*a));
    window.draw(tri);
  }

  int radius = 20;
  sf::CircleShape circ(radius);
  circ.setFillColor(sf::Color::Black);
  circ.setPosition(window.getSize().x/2. - radius, window.getSize().y/2. - radius);
  window.draw(circ);

  window.display();
}

void Visualizer::pollEvents() {
  
  while(window.isOpen()) {
    sf::Event event;

    while( window.pollEvent(event) ) {
      if (event.type == sf::Event::Closed) {
        window.close();
      }
      if (event.type == sf::Event::Resized) {
        window.setView(sf::View(sf::FloatRect(0,0,event.size.width, event.size.height)));
        this -> updateWindow();
      }
    }

  }
}

double Visualizer::getNorm() {
  return norm;
}

void Visualizer::updateNorm(double maxValue) {
  norm = norm + (maxValue - prevMax);
  prevMax = maxValue;
}
