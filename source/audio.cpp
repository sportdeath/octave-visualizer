#include <iostream>
#include <cmath>
#include <complex>

#include <portaudio.h>
#include <fftw3.h>

#include "visualizer.hpp"
#include "audio.hpp"

double AudioStream::quinnKappa(double in) {
    double firstTerm = log(3*in*in + 6*in +1)/4.;
    
    double top = in + 1 - sqrt(2/3.);
    double bot = in + 1 + sqrt(2/3.);
    double secondTerm = sqrt(6)/24. * log(top/bot);
    
    return firstTerm + secondTerm;
}

double AudioStream::quinnsSecondEstimator(int k, 
                                          std::complex<double> left,
                                          std::complex<double> mid,
                                          std::complex<double> right
                                          ) {

    double betaM1 = (left/mid).real();
    double betaP1 = (right/mid).real();

    double deltaM1 = -betaM1/(betaM1 - 1);
    double deltaP1 = betaP1/(betaP1 - 1);
    
    double delta = (deltaM1 + deltaP1)/2. 
                    + quinnKappa(deltaP1 * deltaP1) 
                    - quinnKappa(deltaM1 * deltaM1);
    
    return delta;
}

std::complex<double> fftwToComplex(fftw_complex * fft, int bin) {
  return std::complex<double> (fft[bin][0], fft[bin][1]);
}

void AudioStream::computeSpectrogram() {

  // Clear spectrogram
  for (int i = 0; i < SPEC_SIZE; i++) visualizer -> setSpec(i, 0);

  // Execute FFT
  fftw_execute(fftPlan);

  std::complex<double> left, right, mid;

  // Fill with magnitudes of peak values
  for (int bin = 2; bin < SIZE/2; bin++) {
    left = fftwToComplex(fft, bin - 1);
    mid = fftwToComplex(fft, bin);
    right = fftwToComplex(fft, bin + 1);

    // If the bin is a peak
    if (std::abs(mid) > std::abs(right) && std::abs(mid) > std::abs(left)) {
          
      double delta = quinnsSecondEstimator(bin, left, mid, right);
          
      // If the delta function is a peak
      if ( fabs(delta) < 1) {

        // The peak bin
        double peakBin = bin + delta;
              
        // The peak frequency
        double peakFreq = peakBin * 
                (Pa_GetDeviceInfo(inputParameters.device) -> defaultSampleRate)
                /SIZE;

        // Peak frequency wrapped to an octave between 0 and 1
        // Adjusted so that A = 0 to avoid general splicing between bins
        double AOffset = 1. - fmod(log2(440), 1);
        double freqMod =fmod(log2(peakFreq) + AOffset, 1);
        if (freqMod < 0 ) {
          freqMod += 1;
        }
            
        // The desired position in the octave spectrogram
        double desired = freqMod * SPEC_SIZE;
              
        // Computing the amplitudes across bins and writing to array
        int leftNote = floor(desired);
        int rightNote = (leftNote + 1) % SPEC_SIZE;
              
        double rightPercent = desired - leftNote;
        double leftPercent = 1 - rightPercent;
              
        double rightAmp = rightPercent * std::abs(mid);
        double leftAmp = leftPercent * std::abs(mid);
              
        visualizer -> setSpec(leftNote, visualizer -> getSpec(leftNote) + leftAmp);
        visualizer -> setSpec(rightNote, visualizer -> getSpec(rightNote) + rightAmp);
          }
      }
    }

    //Find max of spectrogram
    double max = visualizer -> getSpec(0);
    for (int i = 0; i < SPEC_SIZE; i++) {
      if ( visualizer -> getSpec(i) > max) {
        max = visualizer -> getSpec(i);
      }
    }

    // Update the normalization factor (HP Filter)
    normalizationFactor += (max - previousMaxValue);
    previousMaxValue = max;

    //Normalize spectrogram using max
    for (int i = 0; i < SPEC_SIZE; i++) {
      if (normalizationFactor != 0) {
        visualizer -> setSpec(i, (visualizer -> getSpec(i))/(normalizationFactor));
      }
    }

    // Update the graphics
    visualizer -> updateWindow();
}

double AudioStream::hammingWindow(int i) {
  return 0.54 - 0.46 * cos( 2 * M_PI * i/((double) SIZE));
}

void AudioStream::processNewSamples ( float * newSamples ) {

  // Shift samples over
  for (int i = 0; i < FRAMES_PER_BUFFER * (OVERLAP - 1); i++) {
    sampleBuffer[i] = sampleBuffer[i + FRAMES_PER_BUFFER];
  }

  // Add new samples
  for (int i = 0; i < FRAMES_PER_BUFFER; i++) {
    sampleBuffer[i + FRAMES_PER_BUFFER * (OVERLAP -1)] = newSamples[i];
  }

  // Window samples
  for (int i = 0; i < SIZE; i++) {
    windowedBuffer[i] = sampleBuffer[i] * hammingWindow(i);
  }

  this -> computeSpectrogram();
}

int AudioStream::visualizerCallback (
    const void * inputBuffer,
    void * outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo * timeInfo,
    PaStreamCallbackFlags statusFlags,
    void * userData ) {

  if ( inputBuffer != NULL ) {

    // Cast user input and data
    AudioStream * s = (AudioStream *) userData;
    float ** channels = ( float ** ) inputBuffer;
    float * input = channels[0];

    // Processing
    s -> processNewSamples(input);
    
  }

  return paContinue;

}

AudioStream::AudioStream() {
  err = Pa_Initialize();
  if( err != paNoError ) this -> streamError();

  int numDevices = Pa_GetDeviceCount();
  if (numDevices ==0) {
    std::cout << "No available audio devices!" << std::endl;
    this -> streamError();
  }

  // List Available devices
  std::cout << "Available devices: " << std::endl;
  for (int i =0; i < numDevices; i++) {
    std::cout << i << ".\t" << 
      Pa_GetDeviceInfo(i) -> name << std::endl;
  }

  // Mark default device
  std::cout << "Device (default = " << 
    Pa_GetDefaultInputDevice() << ". " <<
    Pa_GetDeviceInfo( Pa_GetDefaultInputDevice() ) -> name << 
    "): ";

  // Get user input
  int choice;
  std::string choice_;
  std::getline(std::cin, choice_);
  try {
    choice = std::stoi(choice_);
  } catch (const std::invalid_argument& ia) {
    choice = Pa_GetDefaultInputDevice();
  }

  if (choice < 0 || choice > numDevices) {
    choice = Pa_GetDefaultInputDevice();
  }

  inputParameters.device = choice;

  if (inputParameters.device == paNoDevice) {
    std::cout << "Device " << choice << " not found" << std::endl;
    this -> streamError();
  }

  std::cout << "Using device " << Pa_GetDeviceInfo(choice)-> name << std::endl;
    
  // Initialize input
  inputParameters.channelCount = Pa_GetDeviceInfo(choice) -> maxInputChannels;
  inputParameters.sampleFormat = paFloat32 | paNonInterleaved;
  inputParameters.suggestedLatency = Pa_GetDeviceInfo(choice) -> defaultLowInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;

  // Initialize output
  outputParameters.device = Pa_GetDefaultOutputDevice();
  outputParameters.channelCount = Pa_GetDeviceInfo(outputParameters.device) -> maxOutputChannels;
  outputParameters.sampleFormat = paFloat32;
  outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device ) -> defaultLowOutputLatency;
  outputParameters.hostApiSpecificStreamInfo = NULL;

  // Initialize FFTW
  fftPlan = fftw_plan_dft_r2c_1d(SIZE, 
                              windowedBuffer, 
                              fft,
                              FFTW_ESTIMATE);

  // Initialize normalization
  normalizationFactor = 0;
  previousMaxValue = 0;
}

void AudioStream::startStream( Visualizer * v ) {

  visualizer = v;

  PaError err = Pa_OpenStream (
      &stream,
      &inputParameters,
      &outputParameters,
      Pa_GetDeviceInfo(inputParameters.device) -> defaultSampleRate,
      FRAMES_PER_BUFFER,
      0,
      visualizerCallback,
      this
      );

  if ( err != paNoError ) this -> streamError();

  err = Pa_StartStream( stream );
  if ( err != paNoError ) this -> streamError();
}

void AudioStream::streamError() {
  Pa_Terminate();
  std::cerr << Pa_GetErrorText( err ) << std::endl;
  throw 20;
}

AudioStream::~AudioStream() {
  fftw_destroy_plan(fftPlan);
  Pa_Terminate();
}
