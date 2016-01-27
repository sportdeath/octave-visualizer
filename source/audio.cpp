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

    //double squaredMagnitude = fft[k][0]* fft[k][0] + fft[k][1]* fft[k][1];
    
    //double alpha1 = (fft[k - 1][0]*fft[k][0] + fft[k - 1][1]*fft[k][1])
                    ///squaredMagnitude; //Re(fft[k-1]/fft[k])
    //double alpha2 = (fft[k + 1][0]*fft[k][0] + fft[k + 1][1]*fft[k][1])
                    ///squaredMagnitude; //
    
    double delta1 = alpha1/(1 - alpha1);
    double delta2 = -alpha2/(1 - alpha2);
    
    double delta = (delta1 + delta2)/2. 
                    - quinnKappa(delta1*delta1) 
                    + quinnKappa(delta2 * delta2);
    
    return delta;
}

void computeSpectrogram(double * in, int size, UserData * u) {

  Visualizer * v = u -> visualizer;

  for (int i = 0; i < SPEC_SIZE; i++) v -> setSpec(i, 0);
    
  fftw_complex * fftOut = new fftw_complex[size/2 + 1];
  
  fftw_plan plan = fftw_plan_dft_r2c_1d(size, in, fftOut, FFTW_ESTIMATE);
  
  fftw_execute(plan);
  
  fftw_destroy_plan(plan);

  //Fill with magnitudes of peak values
  for (int bin = 2; bin < size/2; bin++) {
      
      double leftMag = fftOut[bin-1][0]*fftOut[bin-1][0]+fftOut[bin-1][1]*fftOut[bin-1][1];
      double midMag = fftOut[bin][0]*fftOut[bin][0]+fftOut[bin][1]*fftOut[bin][1];
      double rightMag = fftOut[bin+1][0]*fftOut[bin+1][0]+fftOut[bin+1][1]*fftOut[bin+1][1];
      
      if ((midMag >rightMag) && (midMag > leftMag)) {
          
          double delta = quinnsSecondEstimator(bin,fftOut);
          
          if ( fabs(delta) < 1) {
              double peakBin = bin + delta;
              
              double peakFreq = peakBin * 96000/size;
              double freqMod =fmod(log2(peakFreq), 1);
              if (freqMod < 0 ) {
                  freqMod += 1;
              }
              
              double desired =freqMod* SPEC_SIZE;
              
              int leftNote = floor(desired);
              int rightNote = (leftNote + 1) % SPEC_SIZE;
              
              double rightPercent = desired - leftNote;
              double leftPercent = 1 - rightPercent;
              
              double rightAmp = rightPercent * sqrt(fftOut[bin][0] * fftOut[bin][0] + fftOut[bin][1] * fftOut[bin][1]);
              double leftAmp = leftPercent * sqrt(fftOut[bin][0] * fftOut[bin][0] + fftOut[bin][1] * fftOut[bin][1]);
              

              v -> setSpec(leftNote, v -> getSpec(leftNote) + leftAmp);
              v -> setSpec(rightNote, v -> getSpec(rightNote) + rightAmp);
          }
      }
    }

    //Find max of spectrogram
    double max = v -> getSpec(0);
    for (int i = 0; i < SPEC_SIZE; i++) {
      if ( v -> getSpec(i) > max) {
        max = v -> getSpec(i);
      }
    }


    u -> normlizationFactor += (max - u -> prevMaxValue);
    u -> prevMaxValue = max;

    //Normalize spectrogram using max
    for (int i = 0; i < SPEC_SIZE; i++) {
      if (u -> normlizationFactor != 0) {
        v -> setSpec(i, (v -> getSpec(i))/(u -> normlizationFactor));
      }
    }

    v -> updateWindow();
}

double hamming(int i, int N) {
  return 0.54 - 0.46 * cos( 2 * M_PI * i/(double) (N -1));
}

static int visualizerCallBack (
    const void * inputBuffer,
    void * outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo * timeInfo,
    PaStreamCallbackFlags statusFlags,
    void * userData ) {

  if ( inputBuffer != NULL ) {

    float ** channels = ( float ** ) inputBuffer;

    float * in_ = channels[0];

    double in[framesPerBuffer * OVERLAP];

    UserData * u = (UserData *) userData;

    double * prev = u -> previousSamples;

    for (int i = 0; i < FRAMES_PER_BUFFER * (OVERLAP - 1); i++) {
      in[i] = prev[i + FRAMES_PER_BUFFER] * hamming(i, FRAMES_PER_BUFFER * OVERLAP);
      prev[i] = prev[i + FRAMES_PER_BUFFER];
    }

    for (int i = 0; i < framesPerBuffer; i++) {
      int j = i + framesPerBuffer * (OVERLAP -1);
      in[j] = in_[i] * hamming(j, FRAMES_PER_BUFFER * OVERLAP);
      prev[j] = in_[i];
    }

    computeSpectrogram( in, FRAMES_PER_BUFFER * OVERLAP, u );

  }

  return paContinue;

}

portAudioSettings portAudioInit() {
  portAudioSettings s;
  PaError err;

  err = Pa_Initialize();
  if( err != paNoError ) portAudioError(err);

  int numDevices = Pa_GetDeviceCount();

  if (numDevices ==0) portAudioError(err);

  std::cout << "Available devices: " << std::endl;

  for (int i =0; i < numDevices; i++) {
    s.deviceInfo  = Pa_GetDeviceInfo( i );
    std::cout << i << ".\t" << s.deviceInfo -> name << std::endl;
  }

  std::cout << "Device (default = " << 
    Pa_GetDefaultInputDevice() << ". " <<
    Pa_GetDeviceInfo( Pa_GetDefaultInputDevice() ) -> name << 
    "): ";

  int choice;
  std::cin >> choice; 

  if (choice < 0 || choice > numDevices) {
    choice = Pa_GetDefaultInputDevice();
  }

  s.inputParameters.device = choice;
  s.deviceInfo = Pa_GetDeviceInfo( choice );

  std::cout << "Using device " << s.deviceInfo -> name << std::endl;
    
  s.inputParameters.channelCount = s.deviceInfo -> maxInputChannels;
  s.inputParameters.sampleFormat = paFloat32 | paNonInterleaved;
  s.inputParameters.suggestedLatency = s.deviceInfo -> defaultLowInputLatency;
  s.inputParameters.hostApiSpecificStreamInfo = NULL;

  s.outputParameters.device = Pa_GetDefaultOutputDevice();
  s.outputParameters.channelCount = s.deviceInfo -> maxInputChannels;
  s.outputParameters.sampleFormat = paFloat32;
  s.outputParameters.suggestedLatency = Pa_GetDeviceInfo( s.outputParameters.device ) -> defaultLowOutputLatency;
  s.outputParameters.hostApiSpecificStreamInfo = NULL;

  return s;
}

void portAudioStartStream( portAudioSettings s, void * userData ) {

  PaError err = Pa_OpenStream (
      &(s.stream),
      &(s.inputParameters),
      &(s.outputParameters),
      s.deviceInfo -> defaultSampleRate,
      FRAMES_PER_BUFFER,
      0,
      visualizerCallBack,
      userData
      );

  if ( err != paNoError ) portAudioError(err);

  err = Pa_StartStream( s.stream );
  if ( err != paNoError ) portAudioError(err);
}

void portAudioError(PaError err) {
  Pa_Terminate();
  std::cerr << Pa_GetErrorText( err ) << std::endl;
  throw 20;
}
