# Introduction #

The purpose of this page is to present how the signal to noise ratio is computed in the simulator and derive measurements that can help evaluating accurately the merit of such or such signal transmission scheme.

# Additional white Gaussian noise calculation details #

The approach is fairly simple:
  1. The mean signal power is estimated for the whole length of simulation since we produce all signal samples at once
  1. The inverse S/N ratio (the noise to signal ratio) is applied to calculate the noise power
  1. The noise amplitude is calculated as the square root of the noise power
  1. The noise amplitude is further divided by the square root of 2 to account for the complex signal having actually 2 channels (I and Q) so the noise is spanning over twice the space
  1. A unit white Gaussian noise is generated
  1. The unit white Gaussian noise is multiplied by the calculated noise amplitude

# Signal bandwidth considerations #

The simple approach described seems to compute the SNR accurately. However this SNR is taken in the whole complex sampling pass band i.e. the sampling frequency f<sub>s</sub> while we use a coherent correlation scheme with BPSK modulation and this acts as a matching filter ("coherent" is important here) that will mostly ignore the noise outside the signal useful pass band. The useful pass band of the signal which is the actual pass band if the signal is passed through a pulse shaping filter consists of the main lobe around zero frequency from -f<sub>c</sub> to +f<sub>c</sub> so its bandwidth is 2f<sub>c</sub>

Hence the useful noise pass band matches the signal useful pass band and therefore the actual noise power to be considered should be multiplied by 2f<sub>c</sub>/f<sub>s</sub> which means the SNR corresponding to the experimental conditions (for example the limit SNR for near 100% correct decodes) should be multiplied by the inverse that is f<sub>s</sub>/2f<sub>c</sub>

This can be verified experimentally as summarized in the following table. We consider a 1023 bit sequence (m=10) and various chip rates f<sub>c</sub> for various sampling frequencies f<sub>s</sub>:

| **f<sub>s</sub>** | **f<sub>c</sub>** | **f<sub>s</sub>/2f<sub>c</sub>** | **adjust**<br>SNR(dB)<table><thead><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 4096 </td><td> 1023 </td><td> 2 </td><td> 3 </td><td> -24 </td><td> -21 </td></tr>
<tr><td> 8192 </td><td> 1023 </td><td> 4 </td><td> 6 </td><td> -27 </td><td> -21 </td></tr>
<tr><td> 16384 </td><td> 1023 </td><td> 8 </td><td> 9 </td><td> -30 </td><td> -21 </td></tr>
<tr><td> 8192 </td><td> 2046 </td><td> 2 </td><td> 3 </td><td> -24 </td><td> -21 </td></tr>
<tr><td> 4096 </td><td> 511.5 </td><td> 4 </td><td> 6 </td><td> -27 </td><td> -21 </td></tr></tbody></table>

We can see that the SNR per signal bandwidth has a constant value of -21 dB. This appears to be a characteristic of the code.<br>
<br>
For non-coherent correlations with OOK or DBPSK modulations the increase in bandwidth does not yield the corresponding decrease in the limit S/N. For every doubling there is only about 1dB difference. This could be due to the less efficient filtering effect of the correlation.<br>
<br>
<h1>How does this relate to link budget?</h1>

Now we want to calculate the minimum power necessary to reach the limit SNR. We only consider the effects of thermal noise which spectral power is -174 dBm/Hz. Considering the constant SNR per signal bandwidth 2f<sub>c</sub> of -21 dB we have the following for the various chip frequencies:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 511.5 </td><td> 30 </td><td> -144 </td><td> -165 </td></tr>
<tr><td> 1023 </td><td> 33 </td><td> -141 </td><td> -162 </td></tr>
<tr><td> 2046 </td><td> 36 </td><td> -138 </td><td> -159 </td></tr></tbody></table>

Comparatively GPS L1 signals are received on the ground with a power of about -130 dBm for a noise power of -111 dBm. 2f<sub>c</sub> is 2 MHz which represents 63 dB<sub>Hz</sub>: -174 + 63 = -111. Hence the SNR is -19 dB. This  does not exactly match but is still consistent with the -21 dB of the 1023 length code (same as GPS).