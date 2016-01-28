#include <iostream>
#include <cmath>
#include <array>

#include <SFML/Graphics.hpp>

#include "visualizer.hpp"
#include "graphics.hpp"
#include "audio.hpp"

Visualizer::Visualizer(AudioStream * s) {
  stream = s;

  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  settings.antialiasingLevel = 4;
  fullscreen = false;

  window.create(
      sf::VideoMode(600, 600),
      "Visualizer",
      sf::Style::Titlebar |
      sf::Style::Resize |
      sf::Style::Close,
      settings);

  window.setVerticalSyncEnabled(true);
}

void Visualizer::setSpec(int index, double val) {
  spectrogram[index] = val;
}

double Visualizer::getSpec(int index) {
  return spectrogram[index];
}

void Visualizer::updateWindow() {

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
        stream -> ~AudioStream();
        window.close();
      }
      if (event.type == sf::Event::Resized) {
        window.setView(sf::View(sf::FloatRect(0,0,event.size.width, event.size.height)));
        this -> updateWindow();
      }
      if ((event.type == sf::Event::KeyPressed) 
          && sf::Keyboard::isKeyPressed(sf::Keyboard::F)) {
        fullscreen = !fullscreen;
        window.close();
        window.create(
            sf::VideoMode(600, 600),
            "Visualizer",
            fullscreen ?
            sf::Style::Fullscreen :
            sf::Style::Titlebar |
            sf::Style::Resize |
            sf::Style::Close,
            settings);
        window.setVerticalSyncEnabled(true);
      }
    }

  }
}
