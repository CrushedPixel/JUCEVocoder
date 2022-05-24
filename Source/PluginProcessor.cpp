#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
VocoderProcessor::VocoderProcessor()
	: AudioProcessor(BusesProperties()
		.withInput("Input", AudioChannelSet::stereo(), true) // (optional) carrier
        .withInput("Sidechain", AudioChannelSet::stereo(), true) // modulator
		.withOutput("Output", AudioChannelSet::stereo(), true)
	), fftCarrier(fftOrder), fftModulator(fftOrder), fftInverse(fftOrder), receivedMidi(false)
{
    
    mySynth.clearVoices();
    
    for (int i = 0; i < 7; i++)
    {
        mySynth.addVoice (new SynthVoice());
    }
    
    mySynth.clearSounds();
    
    mySynth.addSound (new SynthSound());
    
    
    
}

VocoderProcessor::~VocoderProcessor()
{
}

//==============================================================================
const String VocoderProcessor::getName() const
{
    return JucePlugin_Name;
}

bool VocoderProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool VocoderProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool VocoderProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double VocoderProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int VocoderProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int VocoderProcessor::getCurrentProgram()
{
    return 0;
}

void VocoderProcessor::setCurrentProgram (int index)
{
}

const String VocoderProcessor::getProgramName (int index)
{
    return {};
}

void VocoderProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void VocoderProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    ignoreUnused(samplesPerBlock);
	// Perhaps use this to inform the host of your plug-in's latency
	//setLatencySamples(960);
    lastSampleRate = sampleRate;
    mySynth.setCurrentPlaybackSampleRate(lastSampleRate);

	// Create sine envelop for FFT
	for (int envSamp = 0; envSamp < fftSize; ++envSamp)
		sinenv[envSamp] = sin(((float)envSamp / fftSize) * PI);

	// Zero some buffers
	zeromem(fftCarrierTemp, sizeof(fftCarrierTemp));
	zeromem(fftModulatorTemp, sizeof(fftModulatorTemp));
	zeromem(outStore, sizeof(outStore));

    receivedMidi = false;
}

void VocoderProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool VocoderProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void VocoderProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    if (!midiMessages.isEmpty()) { receivedMidi = true; }

	// TODO: add gates for both MIDI and audio inputs, disallowing any vocoding if either stream is silent
	//		 bandlimit saw synth to get rid of aliasing (will have to use the user-chosen sampling rate to determine cut-off frequency)
	//       somehow obtain lower latency?? Smaller FFT size with better windowing? Different vocoding technique altogether?
	//       replace arrays with std::vector or juce::HeapBlock; This may avoid placing large vectors on the stack, and instead dynamically place them in the heap
	/*
	if (timerFlag) {
		started = std::chrono::high_resolution_clock::now(); // Start timer
		timerFlag = false;
	}
	*/

	ScopedNoDenormals noDenormals;
	auto totalNumInputChannels = getTotalNumInputChannels();
	auto totalNumOutputChannels = getTotalNumOutputChannels();
	auto numSamples = buffer.getNumSamples();

	// In case we have more outputs than inputs, this code clears any output
	// channels that didn't contain input data, (because these aren't
	// guaranteed to be empty - they may contain garbage).
	// This is here to avoid people getting screaming feedback
	// when they first compile a plugin, but obviously you don't need to keep
	// this code if your algorithm always overwrites all the output channels.
	for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
		buffer.clear(i, 0, numSamples);

	synthBuf.clear();
	synthBuf.setSize(totalNumOutputChannels, numSamples);

	mySynth.renderNextBlock(synthBuf, midiMessages, 0, numSamples);
	midiMessages.clear();
    
    // When running in Logic as a MIDI-controlled effect, the sidechain input is provided
    // on bus buffer 0 instead of bus buffer 1. Buffer 1 is empty in that scenario.
    auto carrierBufferIndex = receivedMidi ? 1 : 0;
    auto modulatorBufferIndex = receivedMidi ? 0 : 1;
	
	// Get read pointers to first channel - copy results to both channels (treat as mono)
    auto carrierBuffer = getBusBuffer(buffer, true, carrierBufferIndex);
    auto *carrierChannelData = carrierBuffer.getReadPointer(0);
    
	auto* midiChannelData = synthBuf.getReadPointer(0);
    
    // if there is only a sidechain and no main input, use the sidechain which is then in channel 0
    auto modulatorBuffer = getBusBuffer(buffer, true, modulatorBufferIndex);
    auto *modulatorChannelData = modulatorBuffer.getReadPointer(0);

	for (int sample = 0; sample < numSamples; ++sample)
	{
		carrierInputQueue.push(modulatorChannelData[sample]);
		modulatorInputQueue.push(carrierChannelData[sample] + midiChannelData[sample]);
	}

	// Probably not necessary - but ensures that no original input audio makes it out
	buffer.clear();

	// Compute vocoded result while enough input samples are available
	while (modulatorInputQueue.size() >= fftSize/2)
	{

		// Copy input samples into FFT processing buffer, moving half an FFT block at a time for overlap-add
		for (int i = 0; i < (fftSize / 2); i++)
		{
			fftCarrierTemp[i] = fftCarrierTemp[i + (fftSize / 2)];
			fftModulatorTemp[i] = fftModulatorTemp[i + (fftSize / 2)];

			fftCarrierTemp[i + (fftSize / 2)] = carrierInputQueue.front();
			fftModulatorTemp[i + (fftSize / 2)] = modulatorInputQueue.front();

			if (!carrierInputQueue.empty())
				carrierInputQueue.pop();
			if (!modulatorInputQueue.empty())
				modulatorInputQueue.pop();
		}

		// Change to sizeof(double) if you ever switch to doubles!
		zeromem(fftCarrierOut, sizeof(fftCarrierOut));
		memcpy(fftCarrierOut, fftCarrierTemp, fftSize * sizeof(float));
		zeromem(fftModulatorOut, sizeof(fftModulatorOut));
		memcpy(fftModulatorOut, fftModulatorTemp, fftSize * sizeof(float));

		float inputAudioRMSAmplitude = calculateRMSAmplitudeOfBlock(fftCarrierOut);

		float* vocoderOutput = vocode(fftCarrierOut, fftModulatorOut);

		float outputAudioRMSAmplitude = calculateRMSAmplitudeOfBlock(vocoderOutput);
		
		// prevent division by 0
		if (outputAudioRMSAmplitude > 0)
		{
			// Scale vocoder output using ratio of input and output rms values
			float rmsCorrection = inputAudioRMSAmplitude / outputAudioRMSAmplitude;
			for (int i = 0; i < fftSize; ++i)
			{
				vocoderOutput[i] *= (rmsCorrection * 0.7);

				// There is still occasionally clipping - perhaps add a soft limiter
			}
		}
		
		// Write half block to output queue
		for (int i = 0; i < (fftSize / 2); i++)
		{
			// overlap-add process
			sendOut = (vocoderOutput[i] + outStore[i + (fftSize / 2)]) / 2;
			fftOutputQueue.push(sendOut);
		}

		// store output in buffer for next block (for OLA)
		// Maybe only need to copy half the buffer here - the half that will overlap?
		memcpy(outStore, vocoderOutput, fftSize * sizeof(float));

	}

	//char str[256];
	//sprintf_s(str, "There were %i samples in this block.\n", numSamples);
	//OutputDebugString(str);

	// Write vocoding output to output buffer if enough samples are available
	if (fftOutputQueue.size() >= numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			for (int channel = 0; channel < totalNumInputChannels; ++channel)
			{
				auto* audioChannelDataOut = buffer.getWritePointer(channel);
				audioChannelDataOut[i] = fftOutputQueue.front();
				// For each sample, pop on the last channel written to
				if ((channel == totalNumInputChannels - 1) && !fftOutputQueue.empty())
					fftOutputQueue.pop();
			}
		}
		// Calculate time passed and display in VS17 output
		/*
		if (!timerFlag)
		{
			auto done = std::chrono::high_resolution_clock::now();
			int64_t diffTime = std::chrono::duration_cast<std::chrono::milliseconds>(done - started).count();
			char str[256];
			sprintf_s(str, "This Vocode block took %" PRId64 " milliseconds to calculate.\n", diffTime);
			OutputDebugString(str);
			timerFlag = true;
		}
		*/
	}
		

}

// Make function more general - don't rely on fftSize variable
float VocoderProcessor::calculateRMSAmplitudeOfBlock(float* audioArray)
{
	float sumOfSquares = 0;

	// Use range-based for loop instead
	for (int i = 0; i < fftSize; ++i)
	{
		sumOfSquares += audioArray[i] * audioArray[i];
	}

	float mean = sumOfSquares / fftSize;

	float rms = sqrt(mean);

	return rms;

}

float* VocoderProcessor::vocode(float* modulator, float* carrier)
{
	multiplyBySineEnvelope(modulator);
	multiplyBySineEnvelope(carrier);

	// Perform forward FFTs
	fftCarrier.performRealOnlyForwardTransform(modulator, true);
	fftModulator.performRealOnlyForwardTransform(carrier, true);

	// Get magnitudes of signal spectrum
	getMagnitudeOfInterleavedComplexArray(modulator);

	// Smooth the signal spectrum to get envelope
	smoothSpectrum();

	// Perhaps dynamically allocate instead, or just make a global variable
	float fftOut[fftSize * 2];
	zeromem(fftOut, sizeof(fftOut));

	// Multiply carrier spectrum by signal spectrum envelope
	for (int outSample = 0; outSample < 2 * fftSize; outSample++)
	{
		fftOut[outSample] = carrier[outSample] * fftCarrierEnv[outSample / 2];
	}

	fftInverse.performRealOnlyInverseTransform(fftOut);

	multiplyBySineEnvelope(fftOut);

	return fftOut;
}

void VocoderProcessor::smoothSpectrum()
{
	zeromem(fftCarrierEnv, sizeof(fftCarrierEnv));
	for (int smoothSamp = 0; smoothSamp < fftSize; smoothSamp++) {
		for (int j = 0; j < 10; j++) {
			if ((smoothSamp - j) >= 0) {
				fftCarrierEnv[smoothSamp] += signalMag[smoothSamp - j];
			}
		}
		fftCarrierEnv[smoothSamp] /= 10;
	}
}

void VocoderProcessor::getMagnitudeOfInterleavedComplexArray(float *array)
{
	zeromem(signalMag, sizeof(signalMag));

	for (int i = 0; i < fftSize; ++i)
	{
		signalMag[i] = sqrt(pow(fftCarrierOut[2 * i], 2) + pow(fftCarrierOut[2 * i + 1], 2));
	}
}

void VocoderProcessor::multiplyBySineEnvelope(float *carrier)
{
	for (int i = 0; i < fftSize; ++i)
	{
		carrier[i] *= sinenv[i];
	}
}

//==============================================================================
bool VocoderProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* VocoderProcessor::createEditor()
{
    return new VocoderEditor (*this);
}

//==============================================================================
void VocoderProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void VocoderProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new VocoderProcessor();
}
