#include <iostream>

#include "visualizer.hpp"
#include "audio.hpp"


int main(int argc, char **argv) {

  try {
    AudioStream s;
    Visualizer v(&s);
    s.startStream(&v);
    v.pollEvents();
  } catch (int e) {
    return -1;
  }

  return 0;
}
