# RSSoft vs KV #

This page compares the performance of two soft-decision based Reed-Solomon decoders:
  * The "KV" decoder from Koetter and Vardy's Code Vector Technologies company which is licensed for free for amateur use with WSJT software in a RS(63,12) binary implementation known as KVASD executable.
  * The RSSoft library that is Open Source and hosted also in Google Code as the [rssoft project](https://code.google.com/p/rssoft/) and that will be used with WSGC.

## A picture is worth 1000 words ##

![https://wsgc.googlecode.com/git/img/RSSoft_vs_KV_JT65B.png](https://wsgc.googlecode.com/git/img/RSSoft_vs_KV_JT65B.png)

We see that RSSoft does not perform as well as KV but is off by only 0.5 dB. This is clearly shown by the yellow line that is the red RSSoft's line shifted by -0.5dB which fits nearly perfectly the KV's blue line. In cases where RS(63,12) is not practical KVASD binary becomes unusable and it can be interesting to use RSSoft instead of a hard-decision RS algorithm. Moreover the KVASD binary is not well documented so it is not easy to use. The 0.5 dB may be worth the headache!

The comparison is based on the figure given at the end of the [JT65 description document](http://physics.princeton.edu/pulsar/K1JT/JT65.pdf). To do the best possible comparison we will use WSGC prototype with JT65 source coding and MFSK modulation bypassing the correlation stuff completely and making it very similar to JT65B in WSJT.

WSJT takes 4096 samples at 11025 Hz sampling rate to form a symbol. It takes a 4096 point FFT accordingly to compute the powers in the MFSK frequency bins. A frequency bin is the unit of symbol bandwidth used in WSJT JT65x modes. Each variant uses a multiple of this bandwidth and JT65B takes two of them thus:
  * the symbol bandwidth for JT65B is 2 FFT bins large so (11025/4096)`*`2 = 5.383301 Hz (approximately 5.4 Hz)
  * the symbol period is 4096/11025 = 0.371519 s (approximately 372 ms)

With WSGC MFSK simulation a sampling frequency of 4096 Hz is assumed so with a 4096 point FFT:
  * we can take a symbol bandwidth of 2 Hz and a symbol length of 1s and  have characteristics equivalent to JT65B scaled by 11025/4096 in frequency and the inverse in time. This ratio is about 4.3 dB. To keep things equal 4.3 dB should be added to the SNR of WSGC/MFSK.
  * the noise is calculated as a complex value at 4096 Hz sampling frequency which means that the noise bandwidth is 4096/2500 wider than the 2500 Hz reference bandwidth of WSJT. So another 2.1 dB must be added to WSGC/MFSK SNR figures to make things equal.

In total 6.4 dB are added to the SNR figures of the WSGC prototype simulations before being compared to the JT65B figures.

In a similar manner to the procedure described for JT65B test 1000 simulated transmissions were run for each of the SNR points. The points were shifted by -6.4 dB. For example to calculate the -23 dB point the WSGC simulator is run with a SNR of -29.4 dB. Half dB points were inserted in the steep part of the slope to get a smoother shape.

## Remark on sampling frequency ##

If sampling at 8192 Hz there is no change in the 4.3 dB factor because the symbol timings and bandwidth stay the same. Conversely the noise bandwidth ratio is now 8192/2500 that is twice the one at 4096 Hz sampling rate therefore if the SNR shift is kept at 6.4 dB we should see the 8192 Hz curve (yellow curve) 3 dB lower than the 4096 curve (red curve). This is indeed what is observed:

![https://wsgc.googlecode.com/git/img/MFSK_f_sampling_effect.png](https://wsgc.googlecode.com/git/img/MFSK_f_sampling_effect.png)

Note: if you look at it closely you can notice that the curve decline is sharper at 8192 Hz than at 4096 Hz. This might be due to the difference in the FFT sizes. To accomodate the 2 Hz symbol bandwidth it is 4096 point at 8192 Hz and 2048 point at 4096 Hz.