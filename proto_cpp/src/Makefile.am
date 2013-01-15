bin_PROGRAMS = wsgc_test wsgc_generator

AM_CXXFLAGS = $(FFTW3F_CFLAGS) 

if HAS_CUDA

bin_PROGRAMS += wsgc_cuda

AM_CXXFLAGS += $(CUDA_CFLAGS)

# ls *.cpp -1 | awk '{printf "%s \\\n", $1}'
wsgc_test_SOURCES = CodeModulator_BPSK.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_OOK_detection.cpp \
    ContinuousPhaseCarrier.cpp \
    CorrelationRecord.cpp \
    AutocorrelationRecord.cpp \
    DecisionBox.cpp \
    DecisionBox_Piloted.cpp \
    DecisionBox_Piloted_And_Synced.cpp \
    DecisionBox_Unpiloted.cpp \
    FadingModelClarke.cpp \
    FadingModel.cpp \
    FadingModelNone.cpp \
    FadingModelWatterson.cpp \
    GaussFIR.cpp \
    GoldCodeGenerator.cpp \
    LocalCodes.cpp \
    LocalCodes_Cuda.cu \
    LocalCodes_Host.cpp \
    LocalCodesFFT.cpp \
    LocalCodesFFT_Host.cpp \
    PilotedMessageCorrelator.cpp \
    PilotedMessageCorrelator_Host.cpp \
    PilotedMessageCorrelator_Cuda.cu \
    UnpilotedMessageCorrelator_Host.cpp \
    UnpilotedMultiplePrnCorrelator.cpp \
    Modulation.cpp \
    MultiplePrnCorrelator.cpp \
    MultiplePrnCorrelator_FreqDep.cpp \
    MultiplePrnCorrelator_FreqIndep.cpp \
    PilotedTrainingMessageCorrelator.cpp \
    PilotedTrainingMessageCorrelator_Host.cpp \
    PilotedTrainingMultiplePrnCorrelator.cpp \
    Options.cpp \
    Path.cpp \
    PilotCorrelationAnalyzer.cpp \
    PilotCorrelationRecord.cpp \
    PilotCorrelator.cpp \
    PilotCorrelator_Host.cpp \
    PilotedMultiplePrnCorrelator.cpp \
    SampleSequencer.cpp \
    SimulatedSource.cpp \
    SinglePrnCorrelator_FreqDep.cpp \
    SinglePrnCorrelator_FreqDep_Host.cpp \
    PrnAutocorrelator.cpp \
    PrnAutocorrelator_Host.cpp \
    SourceFFT.cpp \
    SourceFFT_Host.cpp \
    SourceMixer.cpp \
    wsgc_test.cpp \
    WsgcUtils.cpp \
    CudaManager.cu \
    LocalCodesFFT_Cuda.cu \
    PilotCorrelator_Cuda.cu \
    SinglePrnCorrelator_FreqDep_Cuda.cu \
    SourceFFT_Cuda.cu \
    Cuda_Operators.cu

wsgc_test_LDADD = -lrt $(FFTW3F_LIBS) $(NVLDFLAGS)

wsgc_cuda_SOURCES = CudaManager.cu \
	ContinuousPhaseCarrier.cpp \
	SourceFFT.cpp \
	SourceFFT_Cuda.cu \
	SimulatedSource.cpp \
	GoldCodeGenerator.cpp \
	CodeModulator_BPSK.cpp \
	LocalCodesFFT.cpp \
	LocalCodesFFT_Cuda.cu \
	LocalCodes.cpp \
	LocalCodes_Cuda.cu \
    Cuda_Operators.cu \
	test_cuda.cu \
	wsgc_cuda.cpp
	
wsgc_cuda_LDADD = -lrt $(FFTW3F_LIBS) $(NVLDFLAGS)

else

wsgc_test_SOURCES = CodeModulator_BPSK.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_OOK_detection.cpp \
    ContinuousPhaseCarrier.cpp \
    CorrelationRecord.cpp \
    AutocorrelationRecord.cpp \
    DecisionBox.cpp \
    DecisionBox_Piloted.cpp \
    DecisionBox_Piloted_And_Synced.cpp \
    DecisionBox_Unpiloted.cpp \
    FadingModelClarke.cpp \
    FadingModel.cpp \
    FadingModelNone.cpp \
    FadingModelWatterson.cpp \
    GaussFIR.cpp \
    GoldCodeGenerator.cpp \
    LocalCodes.cpp \
    LocalCodes_Host.cpp \
    LocalCodesFFT.cpp \
    LocalCodesFFT_Host.cpp \
    PilotedMessageCorrelator.cpp \
    PilotedMessageCorrelator_Host.cpp \
    UnpilotedMessageCorrelator_Host.cpp \
    UnpilotedMultiplePrnCorrelator.cpp \
    Modulation.cpp \
    MultiplePrnCorrelator.cpp \
    MultiplePrnCorrelator_FreqDep.cpp \
    MultiplePrnCorrelator_FreqIndep.cpp \
    PilotedTrainingMessageCorrelator.cpp \
    PilotedTrainingMessageCorrelator_Host.cpp \
    PilotedTrainingMultiplePrnCorrelator.cpp \
    Options.cpp \
    Path.cpp \
    PilotCorrelationAnalyzer.cpp \
    PilotCorrelationRecord.cpp \
    PilotCorrelator.cpp \
    PilotCorrelator_Host.cpp \
    PilotedMultiplePrnCorrelator.cpp \
    SampleSequencer.cpp \
    SimulatedSource.cpp \
    SinglePrnCorrelator_FreqDep.cpp \
    SinglePrnCorrelator_FreqDep_Host.cpp \
    PrnAutocorrelator.cpp \
    PrnAutocorrelator_Host.cpp \
    SourceFFT.cpp \
    SourceFFT_Host.cpp \
    SourceMixer.cpp \
    wsgc_test.cpp \
    WsgcUtils.cpp

wsgc_test_LDADD = -lrt $(FFTW3F_LIBS)

endif
                    
wsgc_generator_SOURCES = wsgc_generator.cpp \
	CodeModulator_BPSK.cpp \
	ContinuousPhaseCarrier.cpp \
	GoldCodeGenerator.cpp \
    LocalCodesFFT.cpp \
    WsgcUtils.cpp \
    Options.cpp \
    SimulatedSource.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_OOK_detection.cpp \
    CodeModulator_CW_Test.cpp \
    FadingModel.cpp \
    FadingModelNone.cpp \
    FadingModelClarke.cpp \
    FadingModelWatterson.cpp \
    Path.cpp \
    GaussFIR.cpp \
    Modulation.cpp \
    SampleSequencer.cpp \
    SourceMixer.cpp 

wsgc_generator_LDADD = -lrt $(FFTW3F_LIBS)

.cu.o:
	$(NVCC) $(NVCCFLAGS) $(FFTW3F_CFLAGS) -c -o $@ $<