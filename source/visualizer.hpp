#ifndef VISUALIZER
#define VISUALIZER

#include <array>
#include <SFML/Graphics.hpp>

#define SPEC_SIZE 60

class Visualizer {
  private:
    std::array<double, SPEC_SIZE> spectrogram;
    sf::RenderWindow window;

    double norm;
    double prevMax;
  public:
    Visualizer();

    void setSpec(int index, double val);

    double getSpec(int index);

    void updateWindow();

    void pollEvents();

    double getNorm();

    void updateNorm(double maxValue);

};

#endif
