#ifndef AUDIO
#define AUDIO

#define FRAMES_PER_BUFFER 1024
#define OVERLAP 8
#define SIZE (FRAMES_PER_BUFFER * OVERLAP)

#include <portaudio.h>
#include <fftw3.h>


class AudioStream {
  private:
    Visualizer * visualizer;

    PaStreamParameters inputParameters, outputParameters;
    PaStream * stream;
    PaError err;

    double normalizationFactor;
    double previousMaxValue;

    double sampleBuffer[SIZE];
    double windowedBuffer[SIZE];

    fftw_complex fft[SIZE/2 + 1];
    fftw_plan fftPlan;

    // ^^needs destructor!
  
  public:

    AudioStream();

    void startStream(Visualizer * v);

    static int visualizerCallback (
        const void * inputBuffer,
        void * outputBuffer,
        unsigned long framesPerBuffer,
        const PaStreamCallbackTimeInfo * timeInfo,
        PaStreamCallbackFlags statusFlags,
        void * userData );

    void processNewSamples( float * newSamples );

    double hammingWindow(int i);

    void computeSpectrogram();

    double quinnKappa(double in);

    double quinnsSecondEstimator(int k, fftw_complex * fft);

    void streamError();

};

#endif
