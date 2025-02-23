#pragma once
// For debugging in Visual Studio
#include <iostream>
// To recognize "PRId64" for timing; example found here: https://stackoverflow.com/questions/45980608/using-stdchronodurationrep-with-printf-in-32bit-and-64bit-programs
#include <cinttypes>
// For timing
#include <chrono>
// for OutputDebugString()
//#include <windows.h>
#include <stdio.h>
// For sin(), pow(), sqrt() - maybe only need cmath
#include <math.h>
#include <cmath>
#define PI 3.14159265
#include "../JuceLibraryCode/JuceHeader.h"
#include "SynthSound.h"
#include "SynthVoice.h"


//==============================================================================
/**
*/
class VocoderProcessor  : public AudioProcessor
					
{
public:
    //==============================================================================
    VocoderProcessor();
    ~VocoderProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

	float calculateRMSAmplitudeOfBlock(float* audioArray);

	float* vocode(float* modulator, float* carrier);

	void smoothSpectrum();

	void getMagnitudeOfInterleavedComplexArray(float* array);

	void multiplyBySineEnvelope(float* carrier);

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

	enum
	{
		// Change order to lower the latency
		fftOrder = 10,
		fftSize = 1 << fftOrder
	};

private:
	// Synthesizer objects
    Synthesiser mySynth;
    SynthVoice* myVoice;
    
    double lastSampleRate;
	// Stores output of synthesizer
	AudioBuffer<float> synthBuf;
	// The actual FFT objects
	dsp::FFT fftCarrier, fftModulator, fftInverse;
	// Gets passed through FFT and stores output
	float fftCarrierTemp[fftSize], fftModulatorTemp[fftSize];
	float fftCarrierOut[fftSize * 2], fftModulatorOut[fftSize * 2];
	float outStore[fftSize];
	float signalMag[fftSize];
	float fftModulatorEnv[fftSize];
	// Experimenting with queues (adapting deque containers by default)
	std::queue<float> carrierInputQueue, modulatorInputQueue, fftOutputQueue;
	// 1024-point sine envelope
	double sinenv[fftSize];
	// OLA output
	float sendOut;
	//Timer variables
	bool timerFlag = true;
	std::chrono::steady_clock::time_point started;
    
    // Crude way of detecting whether the plugin is running as a MIDI-controlled effect.
    bool receivedMidi;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VocoderProcessor)
};
