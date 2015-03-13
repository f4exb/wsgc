

# Introduction #

This Wiki page describes the two fading models implemented in the C++ prototype. It presents and discuss their advantages and drawbacks. Following is a chapter about how essential fading characteristics are calculated.

# Clarke's fading model #
## Introduction ##
The Clarke's fading model aims at modelling a Rayleigh fading channel with flat fading. While it is considered as computationally inefficient compared to the Jake's simulator its advantage lies in its relative simplicity. Computation time is not a big concern in this simulation since all samples are created at once before the decoding processing starts. The purpose of the model is to implement the fading channel as a FIR filter. For simplicity the excess delays between each cell have a fixed value.
## External links ##
  * [Matlab implementation](http://gaussianwaves.blogspot.com/2011/05/simulation-of-rayleigh-fading-clarkes.html) Mathuranathan Viswanathan Matlab implementation largely reused in a C++ context

# Watterson's fading model #
## Introduction ##
The Watterson fading model is comparatively more elaborate and mainly aims at simulating ionospheric propagation. It does so by simulating magneto-ionic path components as a pair of paths each being modulated by its own independent, complex, bi-variate Gaussian ergodic random process. These are obtained by filtering generated white noise through a Gaussian low pass FIR filter. The width of the filter determines the Doppler spread and also the frequency of fade drops. This model is able to simulate multipath and therefore non flat fading. In addition each path can be frequency shifted which allows the frequency spectrum notches introduced by non-flat fading to vary in time. This more complete model allows to stress test the transmission more accurately and more completely in a larger context than just ionospheric propagation by choosing path variables appropriately. Therefore this model will be preferred in the tests with fading.
## External links ##
  * [PathSim technical details](http://www.moetronix.com/ae4jy/files/pathsimtech100.pdf) Moe Wheatley, AE4JY _Path Sim_ simulation software details. The fading engine is essentially built using the same technique and large parts of the code have been reused in a C++ context
  * [May-June 2000 QEX Magazine](http://om6bb.bab.sk/files/HAM%20kniznica/Magaziny/QEX/03%20May-June%202000%20QEX.pdf) See Johan B. Forrer, KC7WW article "A Low-Cost HF Channel Simulator for Digital Systems"

# Fading simulation method #

**Process without channel fading (LOS) :**
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_no_fading_process.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_no_fading_process.png)


**Process with channel fading :**
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_fading_process.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_fading_process.png)

There are two main blocks in the fading simulation. The fading model properly speaking as described in the previous two chapters followed by an Additive White Gaussian Noise (AWGN) stage. AWGN can be used exclusively to simulate channels in LOS situations such as direct space communications or short range unobstructed terrestrial communications with no multipath effects. That is an ideal case and of course yields best results however this is rarely a concrete case in terrestrial  communications and never the case for EME.

## AWGN calculation ##
The amount of AWGN is calculated by calculating the signal mean power first and adjusting noise amplitude accordingly. As the complete signal is generated first this can be calculated over the whole simulation time.

Noise amplitude is calculated as the square root of mean signal power divided by linear SNR then divided by the square root of 2 because it is a complex signal and thus affecting both the in phase and quadrature channels.

Noise is calculated as a unity normal distribution then multiplied by noise amplitude. For details check [FadingModel::apply\_awgn](http://f4exb.free.fr/wsgc/proto_cpp/doc/html/classFadingModel.html#a0905c6bea8794dcaf2e16f7d4c487c5b) method in the prototype.

If any channel fading is active the mean signal power is calculated on the faded signal i.e. after applying the Clarke or Watterson fading methods.

AWGN is sampled at the main sampling frequency f<sub>s</sub> and therefore corresponds to the pass-band f<sub>s</sub>/2. However since no filtering is applied consecutively this essentially corresponds to the noise in the pass-band of the signal if the signal was previously filtered to match its main lobe which is what would be done in a real life system.

## Typical fading characteristics for VHF and up ##
Essential fading characteristics are given by the time and frequency spreading values. For flat fading the time spreading is always zero and the Clarke's model can be applied. Conversely the Watterson's model allows the tuning of both parameters. The following table is extracted from Klaus Von Der Heide, DJ5HG's article on weak signal communications published in Dubus magazine 1/2009:

| **Channel** | **Time spreading** | **Frequency spreading** | **Comments** |
|:------------|:-------------------|:------------------------|:-------------|
| Tropo inversion | very low | very low | Almost just AWGN: best case |
| Tropo scatter | 10 us | 5 Hz | Need a fairly fast chip rate. Chip rate is not too much limited by time spreading |
| Aurora | 1 ms | 1 kHz | Good luck! Time and frequency spreading are too large to make any compromise. Try Morse Code!|
| EME 144 MHz | 0.1 ms | 0.2 Hz | Not so bad indeed. Gets worse rapidly as frequency increases (2) |
| Airplane scatter (1) | 10 us | 50 Hz | Need a very fast chip rate. There is also the problem of varying Doppler shift and short burst times |
| ISS scatter | 0 | very low | Also an ideal case. Mind the Doppler shift that can be calculated however |

(1): This data is dubious. Why would it be so much different from ISS apart from the fact that the Doppler shift is difficult to predict? On the few cases I witnessed frequency spreading was not larger than with tropo.

(2): F2CT, Guy reports 200 Hz Doppler spread at 24 GHz. That doesn't make it any worse than optimum Rain Scatter at that frequency. However the delay spread must be much larger than with RS.

Not explicitly mentioned in this table is the case of Rain Scatter. Although CW signals may sound like aurora and do have a large frequency spreading the time spreading is much smaller. This is because the various paths go through the same rain cell and cannot have a large delay with one another particularly if the stations are aligned with the rain cell (best conditions). So basically with a very fast chip rate the message would be able to pass through. The Doppler shift can be fairly large and difficult to predict. It is introduced by the movement of the rain cell over ground. If its path is perpendicular to the direction between stations it will be minimal (again these are the best conditions).

With the Watterson's model path delay hence time spreading and Doppler spread hence frequency spreading can be adjusted.

## Typical fading characteristics at HF ##
The best suitable model is the Watterson's model which is explicitly dedicated to ionospheric propagation. The following table is compiled from ITU-R F.1487 Channel parameters:
| **Latitude** | **Condition** | **Delay (ms)** | **Doppler spread (Hz)** |
|:-------------|:--------------|:---------------|:------------------------|
| Low      | Quiet | 0.5 | 0.5 |
| Low      | Moderate | 2 | 1.5 |
| Low      | Disturbed | 6 | 10 |
| Mid      | Quiet | 0.5 | 0.1 |
| Mid      | Moderate | 1 | 0.5 |
| Mid      | Disturbed | 2 | 1 |
| High     | Quiet | 1 | 0.5 |
| High     | Moderate | 3 | 10 |
| High     | Disturbed | 7 | 30 |

Except in the case of low and high latitude disturbed conditions the essential parameter to mitigate appears to be the time spreading.

