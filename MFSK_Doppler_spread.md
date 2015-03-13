# Doppler spread study #



The purpose is to study the performance of the various codes primarily based on MFSK but we will also introduce briefly GC correlated schemes.

In all the test we do not take into account multipath condition. In fact this is negligible for any practical situation except maybe aurora since the symbol length is at least 0.186s on JT scale (0.5s for WSGC prototype).

We will be using the Watterson model with only one path. This fading model multiplies samples with a random complex value hence inducing both phase and amplitude fading. This random variable follows a Gaussian distribution which sigma is the Doppler spread frequency.

## Small Doppler spread ##

In this section we will deal with a maximum Doppler spread of 2 Hz which in JT65 terms is already 5.4 Hz typical of troposcatter at VHF.

In all figures the line style relates to the 3 Doppler spread values used:
  * Continuous: without Doppler fading
  * Fine dashed: with 1Hz Doppler fading
  * Larger dashes: with 2 Hz Doppler fading

### MFSK modes (JT65 clones) ###

Here we deal with MFSK modes using RS(63,12) that are clones of the JT65 modes. We have the following equivalence taking 6.4 dB difference in SNR because of time and bandwidth differences. See [RSSoft vs KV](https://code.google.com/p/wsgc/wiki/RSSoft_vs_KV) for details. We have also a 0.5 dB performance difference between RSSoft and KVASD decoders:
| JT65 | Symbol<br>BW (Hz) <table><thead><th> Symbol<br> time (s) </th><th> MFSK </th><th> Symbol<br>BW (Hz) </th><th> Symbol<br> time (s) </th></thead><tbody>
<tr><td> JT65A </td><td> 2.7 </td><td> 0.372 </td><td> MFSK:0,0 </td><td> 1 </td><td> 1 </td></tr>
<tr><td> JT65B </td><td> 5.4 </td><td> 0.372 </td><td> MFSK:1,0 </td><td> 2 </td><td> 1 </td></tr>
<tr><td> JT65C </td><td> 10.8 </td><td> 0.372 </td><td> MFSK:2,0 </td><td> 4 </td><td> 1 </td></tr>
<tr><td> N/A </td><td> 21.6 </td><td> 0.372 </td><td> MFSK:3,0 </td><td> 8 </td><td> 1 </td></tr></tbody></table>

<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_0-2_MFSK.png' />

We have the following color codes:<br>
<ul><li><font color='#B00000'><b>Red: MFSK:1,0</b></font>: 2 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65B with 5,4 Hz bandwidth/<br>
</li><li><font color='#00A000'><b>Green: MFSK:2,0</b></font>: 4 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65C with 10.8 Hz bandwidth.<br>
</li><li><font color='#F7ED01'><b>Yellow: MFSK:3,0</b></font>: 8 Hz symbol bandwidth and 1s symbol time. Has no JT65 equivalent, this would be a JT65D with 21.6 Hz bandwidth.</li></ul>

Not surprisingly a larger bin size can better accomodate a larger Doppler spread. For a 2 Hz Doppler spread, degradation is more than 5 dB for a 2 Hz bin (red lines), 2 dB for a 4 Hz bin (green lines) and 0.5 dB for a 8 Hz bin (yellow lines). The price to pay is a performance hit at lower Doppler spreads.<br>
<br>
MFSK:1,0 and MFSK:2,0 at 1 Hz on one hand (the red and green fine dashed lines) and MFSK:2,0 and MFSK:3,0 at 2Hz on the other hand (the gren and yellow larger dashed lines) look very close with each other. In each case you have one mode at a Doppler spread of half its bin size and the other at a Doppler spread of a fourth its bin size but with twice the bin size. In each of the groups it seems that the smaller bin size wins (so with a Doppler spread of half the bin size). This might have to be checked in each case however since the lines are very close but it appears these are the two best options.<br>
<br>
<h3>JT33 modes</h3>

So called JT33 modes are MFSK with a RS(31,15) code. The numbers after the colon relate to their MFSK equivalent.<br>
<br>
<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_0-2_JT33.png' />

Colors:<br>
<ul><li><font color='#B00000'><b>Red: MFSK:1,0</b></font>: 2 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65B with 5,4 Hz bandwidth/<br>
</li><li><font color='#FFB1C9'><b>Salmon: MFSK:2,0</b></font>: 4 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65C with 10.8 Hz bandwidth.<br>
</li><li><font color='#008000'><b>Dark Green: JT33:1,1</b></font>: Uses a RS(31,15) code with 2 Hz symbol bandwidth and 2s symbol time.<br>
</li><li><font color='#B5D46E'><b>Light Green: JT33:2,1</b></font>: Uses a RS(31,15) code with 4 Hz symbol bandwidth and 2s symbol time.</li></ul>

JT33 appears to be affected by Doppler spread as can be expected but not much worse than its MFSK JT65 like counterparts if the Doppler spread is limited to half the symbol bandwidth. Since it does better at low on no Doppler spread it can still be a good choice.<br>
<br>
<h3>JT129 modes</h3>

JT129 modes pack the 72 bits of the original JT source coded message along with 10 more bits for the two last digits of a 6 character locator and 2 spare bits. This makes 84 bits in total that is 12 symbols of 7 bits. Therefore a RS(127,12) code is used and the codeword is 127 symbols long.<br>
<br>
<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_0-2_JT129.png' />

Colors:<br>
<ul><li>Red: MFSK:1,0: 2 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65B with 5,4 Hz bandwidth/<br>
</li><li>Salmon: MFSK:2,0: 4 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65C with 10.8 Hz bandwidth.<br>
</li><li>Dark Green: JT129:1,0: Uses a RS(127,12) code with 2 Hz symbol bandwidth, 1s symbol time and 127 symbols transmitted.<br>
</li><li>Light Green: JT129:2,0: Uses a RS(127,12) code with 4 Hz symbol bandwidth, 1s symbol time and 127 symbols transmitted.<br>
</li><li>Cyan: JT129:2,-1: Uses a RS(127,12) code with 4 Hz symbol bandwidth, 0.5s symbol time and 127 symbols transmitted resulting in 24s message time on JT scale (doubled to 48s with a synchro). This is an attempt at reducing message length and have some Doppler immunity with the 4Hz symbol bandwidthh.</li></ul>

<h3>JT257 modes</h3>

JT257 packs the 72 bits of the original JT encoded source message into 9 symbols of 8 bits with no spare. Therefore it uses a RS(255,9) code with a 255 symbols codeword.<br>
<br>
<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_0-2_JT257.png' />

Colors:<br>
<ul><li>Red: MFSK:1,0: 2 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65B with 5,4 Hz bandwidth/<br>
</li><li>Salmon: MFSK:2,0: 4 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65C with 10.8 Hz bandwidth.<br>
</li><li>Dark Green: JT257:1,0: Uses a RS(255,9) code with 2 Hz symbol bandwidth, 1s symbol time and 255 symbols transmitted.<br>
</li><li>Light Green: JT257:2,0: Uses a RS(255,9) code with 4 Hz symbol bandwidth, 1s symbol time and 255 symbols transmitted.</li></ul>

With 4 Hz symbol bandwidth 1 Hz Doppler yields better results than no Doppler. This is also observed on the JT129 figure. The only explanation I can think of is that the dithering introduced by Doppler spread improves the filling of the frequency bin. This could relate to the bin window shape so the appropriate choice of the FFT window seems to matter at least for JT129 and JT257 modes.<br>
<br>
<h3>GC correlation</h3>

<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_0-2_GC.png' />

Colors:<br>
<ul><li>Red: MFSK:1,0: 2 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65B with 5,4 Hz bandwidth/<br>
</li><li>Salmon: MFSK:2,0: 4 Hz symbol bandwidth and 1s symbol time. This is the equivalent of JT65C with 10.8 Hz bandwidth.<br>
</li><li>Dark green: GC:1023,1,4: BPSK GC correlation with 1023 Hz chip frequency, 1s PRN time (therefore using 1023 chips) and 4 PRNs per symbol.<br>
</li><li>Light green: GC:1023,0.5,4: BPSK GC correlation with 2046 Hz chip frequency, 0.5s PRN time (using 1023 chips) and 4 PRNs per symbol.</li></ul>

As expected it is the worst performer there. The smallest Doppler spread makes it go haywire. An attempt to use a shorter time for the PRN makes it just worse even with no Doppler. A pity this was the main subject of this work but it proves again it is completely unusable for weak signals transmission.<br>
<br>
<h2>Large Doppler spread</h2>

We will study the following symbol bandwidths mostly with 1s symbol time (n=0):<br>
<br>
<table><thead><th> BW (Hz) </th><th> Mode </th><th> JT BW (Hz) </th><th> 1/2 JT BW (Hz) </th></thead><tbody>
<tr><td> 32 </td><td> :5,n </td><td> 86 </td><td> 43 </td></tr>
<tr><td> 64 </td><td> :6,n </td><td> 172 </td><td> 86 </td></tr>
<tr><td> 128 </td><td> :7,n </td><td> 344 </td><td> 172 </td></tr>
<tr><td> 256 </td><td> :8,n </td><td> 688 </td><td> 344 </td></tr></tbody></table>

Bandwidths of 344 and 172 Hz in JT terms are close to the JT4G and JT4F bandwidths. For low Doppler spread we have seen that we can expect a  severe drop in performance when the Doppler spread is larger than half of the symbol bandwidth. However with large Doppler it seems to be different so we will study half and full symbol bandwidth Doppler spreads.<br>
<br>
In all figures the line style relates to the 2 Doppler spread values used:<br>
<ul><li>Continuous: without Doppler fading<br>
</li><li>Fine dashed: with half symbol bandwidth Doppler fading<br>
</li><li>Larger dash: with symbol bandwidth Doppler fading</li></ul>

<h3>MFSK</h3>

<img src='https://wsgc.googlecode.com/git/img/MFSK_DopplerSpread_large_MFSK.png' />

Colors:<br>
<ul><li>Blue: JT65B reference<br>
</li><li>Red: MFSK:5,0: 32 Hz symbol bandwidth and 1s symbol time.<br>
</li><li>Green: MFSK:6,0: 64 Hz symbol bandwidth and 1s symbol time.<br>
</li><li>Cyan: MFSK:7,0: 128 Hz symbol bandwidth and 1s symbol time.<br>
</li><li>Deep violet: MFSK:8,0: 256 Hz symbol bandwidth and 1s symbol time.