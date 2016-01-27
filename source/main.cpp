#include <iostream>

#include "visualizer.hpp"
#include "audio.hpp"

// To Do
// Clean up
// exceptions
// Windowing


int main(int argc, char **argv) {

  portAudioSettings s;
  try {
    s = portAudioInit();
  } catch (int e) {
    return -1;
  }

  UserData u;
  Visualizer v;
  u.visualizer = &v;
  u.normlizationFactor = 0;
  u.prevMaxValue = 0;
  for (int i = 0; i < FRAMES_PER_BUFFER * OVERLAP; i++) {
    u.previousSamples[i] = 0;
  }

  try {
    portAudioStartStream(s, &u);
  } catch (int e) {
    return -1;
  }

  v.updateWindow();

  v.pollEvents();

  return 0;
}
