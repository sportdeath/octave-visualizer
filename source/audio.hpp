#ifndef AUDIO
#define AUDIO

#define FRAMES_PER_BUFFER 1024
#define OVERLAP 8

#include <portaudio.h>
#include <fftw3.h>


double quinnKappa(double in);

double quinnsSecondEstimator(int k, fftw_complex * fft);

void computeSpectrogram(double * in, int size, Visualizer * v);

struct portAudioSettings {
  PaStreamParameters inputParameters, outputParameters;
  PaStream * stream;
  const PaDeviceInfo * deviceInfo;
};

struct UserData {
  Visualizer * visualizer;
  double normlizationFactor;
  double prevMaxValue;
  double previousSamples[FRAMES_PER_BUFFER * OVERLAP];
};

portAudioSettings portAudioInit();

void portAudioStartStream(portAudioSettings s, void * userData);

void portAudioError( PaError err );

#endif
