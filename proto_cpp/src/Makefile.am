bin_PROGRAMS = wsgc_test wsgc_generator wsgc_jt65src
wsgc_test_SOURCES = Options.cpp
wsgc_generator_SOURCES = Options.cpp
wsgc_test_LDADD =
wsgc_generator_LDADD =

AM_CXXFLAGS = $(FFTW3F_CFLAGS) $(BOOST_CPPFLAGS)
RSSoft_Engine.$(OBJEXT) : CXXFLAGS += -std=c++0x 

if HAS_RSSOFT
wsgc_test_SOURCES += RSSoft_Engine.cpp \
    RSSoft_DecisionBox.cpp
wsgc_test_LDADD += $(RSSOFT_LDFLAGS)
wsgc_generator_SOURCES += RSSoft_Engine.cpp
wsgc_generator_LDADD += $(RSSOFT_LDFLAGS)
AM_CXXFLAGS += $(RSSOFT_CFLAGS)
NVCCFLAGS += $(RSSOFT_CFLAGS)
endif

if HAS_CUDA

bin_PROGRAMS += wsgc_cuda

AM_CXXFLAGS += $(CUDA_CFLAGS)

# ls *.cpp -1 | awk '{printf "%s \\\n", $1}'
wsgc_test_SOURCES += CodeModulator_BPSK.cpp \
    CodeModulator_DBPSK.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_MFSK.cpp \
    ContinuousPhaseCarrier.cpp \
    TimeCorrelationAnalyzer.cpp \
    CorrelationRecord.cpp \
    AutocorrelationRecord.cpp \
    DecisionRecord.cpp \
    DecisionBox.cpp \
    DecisionBox_Piloted_And_Synced.cpp \
    DecisionBox_Unpiloted_And_Synced.cpp \
    DecisionBox_Thresholds.cpp \
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
    Modulation.cpp \
    DemodulatorSquaring.cpp \
    DemodulatorDifferential.cpp \
    PilotedTrainingMessageCorrelator.cpp \
    PilotedTrainingMessageCorrelator_Host.cpp \
    PilotedTrainingMessageCorrelator_Cuda.cu \
    PilotedTrainingMultiplePrnCorrelator.cpp \
    TrainingCorrelationRecord.cpp \
    MFSK_Options.cpp \
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
    UnpilotedMultiplePrnCorrelator.cpp \
    UnpilotedMultiplePrnCorrelator_Host.cpp \
    UnpilotedTrainingMessageCorrelator.cpp \
    UnpilotedTrainingMessageCorrelator_Host.cpp \
    SourceFFT.cpp \
    SourceFFT_Host.cpp \
    SourceMixer.cpp \
    wsgc_test.cpp \
    WsgcUtils.cpp \
    FIR_RCoef.cpp \
    FIRCoefGenerator.cpp \
    FIRCoefGenerator_RCos.cpp \
    MFSK_MessageDemodulationRecord.cpp \
    MFSK_MessageDemodulator.cpp \
    MFSK_MessageDemodulator_Host.cpp \
    DecisionBox_MFSK.cpp \
    SourceCodec_JT65.cpp \
    Locator.cpp \
    CudaManager.cu \
    CudaDeviceManager.cu \
    LocalCodesFFT_Cuda.cu \
    PilotCorrelator_Cuda.cu \
    SinglePrnCorrelator_FreqDep_Cuda.cu \
    SourceFFT_Cuda.cu \
    Cuda_Operators.cu 
if HAS_RSSOFT
wsgc_test_SOURCES += RSSoft_ReliabilityMatrixCuda.cu
endif	
    
wsgc_cuda_SOURCES = CudaManager.cu \
	CudaDeviceManager.cu \
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
	
wsgc_test_LDADD += $(NVLDFLAGS)
wsgc_cuda_LDADD = -lrt $(FFTW3F_LIBS) $(NVLDFLAGS)	
	
else

wsgc_test_SOURCES += CodeModulator_BPSK.cpp \
    CodeModulator_DBPSK.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_MFSK.cpp \
    ContinuousPhaseCarrier.cpp \
    TimeCorrelationAnalyzer.cpp \
    CorrelationRecord.cpp \
    AutocorrelationRecord.cpp \
    DecisionRecord.cpp \
    DecisionBox.cpp \
    DecisionBox_Piloted_And_Synced.cpp \
    DecisionBox_Unpiloted_And_Synced.cpp \
    DecisionBox_Thresholds.cpp \
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
    Modulation.cpp \
    DemodulatorSquaring.cpp \
    DemodulatorDifferential.cpp \
    PilotedTrainingMessageCorrelator.cpp \
    PilotedTrainingMessageCorrelator_Host.cpp \
    PilotedTrainingMultiplePrnCorrelator.cpp \
    TrainingCorrelationRecord.cpp \
    MFSK_Options.cpp \
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
    UnpilotedMultiplePrnCorrelator.cpp \
    UnpilotedMultiplePrnCorrelator_Host.cpp \
    UnpilotedTrainingMessageCorrelator.cpp \
    UnpilotedTrainingMessageCorrelator_Host.cpp \
    SourceFFT.cpp \
    SourceFFT_Host.cpp \
    SourceMixer.cpp \
    wsgc_test.cpp \
    WsgcUtils.cpp \
    FIR_RCoef.cpp \
    FIRCoefGenerator.cpp \
    FIRCoefGenerator_RCos.cpp \
    MFSK_MessageDemodulationRecord.cpp \
    MFSK_MessageDemodulator.cpp \
    MFSK_MessageDemodulator_Host.cpp \
    DecisionBox_MFSK.cpp \
    SourceCodec_JT65.cpp \
    Locator.cpp

endif

wsgc_generator_SOURCES += wsgc_generator.cpp \
	CodeModulator_BPSK.cpp \
	CodeModulator_DBPSK.cpp \
	ContinuousPhaseCarrier.cpp \
	GoldCodeGenerator.cpp \
    LocalCodesFFT.cpp \
    WsgcUtils.cpp \
    MFSK_Options.cpp \
    DecisionBox_Thresholds.cpp \
    SimulatedSource.cpp \
    CodeModulator_OOK.cpp \
    CodeModulator_CW_Test.cpp \
    FadingModel.cpp \
    FadingModelNone.cpp \
    FadingModelClarke.cpp \
    FadingModelWatterson.cpp \
    Path.cpp \
    GaussFIR.cpp \
    Modulation.cpp \
    SampleSequencer.cpp \
    SourceMixer.cpp \
    FIR_RCoef.cpp \
    FIRCoefGenerator.cpp \
    FIRCoefGenerator_RCos.cpp \
    DemodulatorDifferential.cpp \
    DemodulatorSquaring.cpp \
    CodeModulator_MFSK.cpp \
    SourceCodec_JT65.cpp \
    Locator.cpp
    
wsgc_jt65src_SOURCES = wsgc_jt65src.cpp \
    SourceCodec_JT65.cpp \
    Locator.cpp
    
wsgc_jt65src_LDADD = -lrt -lm
wsgc_generator_LDADD += -lrt $(FFTW3F_LIBS)
wsgc_test_LDADD += -lrt $(FFTW3F_LIBS)

.cu.o:
	$(NVCC) $(NVCCFLAGS) $(FFTW3F_CFLAGS) -c -o $@ $<