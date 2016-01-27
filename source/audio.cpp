#include <iostream>
#include <cmath>
#include <complex>

#include <portaudio.h>
#include <fftw3.h>

#include "visualizer.hpp"
#include "audio.hpp"

double quinnKappa(double in) {
    double firstTerm = log(3*in*in + 6*in +1)/4.;
    
    double top = in + 1 - sqrt(2/3.);
    double bot = in + 1 + sqrt(2/3.);
    double secondTerm = sqrt(6)/24. * log(top/bot);
    
    return firstTerm + secondTerm;
}

double quinnsSecondEstimator(int k, fftw_complex * fft) {

    std::complex<double> present(fft[k][0], fft[k][1]);
    std::complex<double> past(fft[k - 1][0], fft[k - 1][1]);
    std::complex<double> future(fft[k + 1][0], fft[k + 1][1]);
    
    double alpha1 = (past/present).real();
    double alpha2 = (future/present).real();

    double delta1 = alpha1/(1 - alpha1);
    double delta2 = -alpha2/(1 - alpha2);
    
    double delta = (delta1 + delta2)/2. 
                    - quinnKappa(delta1*delta1) 
                    + quinnKappa(delta2 * delta2);
    
    return delta;
}

void AudioStream::computeSpectrogram() {

  for (int i = 0; i < SPEC_SIZE; i++) visualizer -> setSpec(i, 0);
  
  fftw_execute(fftPlan);

  //Fill with magnitudes of peak values
  for (int bin = 2; bin < SIZE/2; bin++) {
      
      double leftMag = fft[bin-1][0]*fft[bin-1][0]+fft[bin-1][1]*fft[bin-1][1];
      double midMag = fft[bin][0]*fft[bin][0]+fft[bin][1]*fft[bin][1];
      double rightMag = fft[bin+1][0]*fft[bin+1][0]+fft[bin+1][1]*fft[bin+1][1];
      
      if ((midMag >rightMag) && (midMag > leftMag)) {
          
          double delta = quinnsSecondEstimator(bin,fft);
          
          if ( fabs(delta) < 1) {
              double peakBin = bin + delta;
              
              double peakFreq = peakBin * 
                (Pa_GetDeviceInfo(inputParameters.device) -> defaultSampleRate)
                /SIZE;

              double freqMod =fmod(log2(peakFreq), 1);
              if (freqMod < 0 ) {
                  freqMod += 1;
              }
              
              double desired =freqMod* SPEC_SIZE;
              
              int leftNote = floor(desired);
              int rightNote = (leftNote + 1) % SPEC_SIZE;
              
              double rightPercent = desired - leftNote;
              double leftPercent = 1 - rightPercent;
              
              double rightAmp = rightPercent * sqrt(fft[bin][0] * fft[bin][0] + fft[bin][1] * fft[bin][1]);
              double leftAmp = leftPercent * sqrt(fft[bin][0] * fft[bin][0] + fft[bin][1] * fft[bin][1]);
              

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


    normalizationFactor += (max - previousMaxValue);
    previousMaxValue = max;

    //Normalize spectrogram using max
    for (int i = 0; i < SPEC_SIZE; i++) {
      if (normalizationFactor != 0) {
        visualizer -> setSpec(i, (visualizer -> getSpec(i))/(normalizationFactor));
      }
    }

    visualizer -> updateWindow();
}

double AudioStream::hammingWindow(int i) {
  return 0.54 - 0.46 * cos( 2 * M_PI * i/(double) (FRAMES_PER_BUFFER * OVERLAP));
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
  for (int i = 0; i < FRAMES_PER_BUFFER * OVERLAP; i++) {
    windowedBuffer[i] = hammingWindow(sampleBuffer[i]);
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

    AudioStream * s = (AudioStream *) userData;

    float ** channels = ( float ** ) inputBuffer;

    float * input = channels[0];
    
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

  int choice;
  std::string choice_;
  std::getline(std::cin, choice_);
  try {
    choice = std::stoi(choice_);
  } catch (int e) {
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
    
  inputParameters.channelCount = Pa_GetDeviceInfo(choice) -> maxInputChannels;
  inputParameters.sampleFormat = paFloat32 | paNonInterleaved;
  inputParameters.suggestedLatency = Pa_GetDeviceInfo(choice) -> defaultLowInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;

  outputParameters.device = Pa_GetDefaultOutputDevice();
  outputParameters.channelCount = Pa_GetDeviceInfo(outputParameters.device) -> maxOutputChannels;
  outputParameters.sampleFormat = paFloat32;
  outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device ) -> defaultLowOutputLatency;
  outputParameters.hostApiSpecificStreamInfo = NULL;

  fftw_plan plan = fftw_plan_dft_r2c_1d(SIZE, 
                                        windowedBuffer, 
                                        fft,
                                        FFTW_ESTIMATE);
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
