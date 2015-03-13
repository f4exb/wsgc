**Spectrum and spectrum display considerations**



# Introduction #

The purpose of this page is to show an example of a WSGC signal spectrum and what it looks like when buried in noise.

# GNU Radio graph #

_[GNU Radio 3.6.0](http://gnuradio.org/redmine/news/13) is used here_

WSGC generator is used to generate a baseband I/Q file then this file can be displayed on a spectrum and waterfall QT display in GNU Radio using the following graph:

![https://wsgc.googlecode.com/git/img/wsgc_GnuRadio_graph.png](https://wsgc.googlecode.com/git/img/wsgc_GnuRadio_graph.png)

_FFT default size is 2048 but it can be adjusted in the QT widget_

# Signal spectrum #

## 1023 Hz chip frequency ##

This is a standard BPSK spectrum pulse shaped by a raised cosine FIR filter. Filter has 13 taps as specified in the -L parameter. Here we simulate a training sequence starting at PRN #3 for 40 PRNs. Signal is centered at 2.12 Hz but this is hardly noticeable (this is to exercise the frequency tracking!). Command line is:

`wsgc_generator -s 4096 -c 1023 -C 1020 -t 2.12 -N 4 -I 3 -R 40 -P 2 -L "rcos:13,1.0,2048.0" -o test_BPSK.raw --simulate-trn`

**The spectrum looks like:**

![https://wsgc.googlecode.com/git/img/wsgc_4096_1023_noawgn_nofad_trn.spec.png](https://wsgc.googlecode.com/git/img/wsgc_4096_1023_noawgn_nofad_trn.spec.png)

**And the corresponding waterfall:**

![https://wsgc.googlecode.com/git/img/wsgc_4096_1023_noawgn_nofad_trn.wf.png](https://wsgc.googlecode.com/git/img/wsgc_4096_1023_noawgn_nofad_trn.wf.png)

## 63 Hz chip frequency ##

We will also try with a much smaller bandwidth signal in the hope it will better show in the noise. We reduce the sampling frequency to 1024 Hz to keep FFT in the same range.

`wsgc_generator -s 1024 -c 63 -m 6 -M 48 -C 50 -t 2.12 -N 4 -I 1 -R 36 -P 2 -L "rcos:21,1.0,128.0" -o test_BPSK.raw --simulate-trn`

**The spectrum looks like:**

![https://wsgc.googlecode.com/git/img/wsgc_1024_63_noawgn_nofad_trn.spec.png](https://wsgc.googlecode.com/git/img/wsgc_1024_63_noawgn_nofad_trn.spec.png)

**And the corresponding waterfall:**

![https://wsgc.googlecode.com/git/img/wsgc_1024_63_noawgn_nofad_trn.4k.png](https://wsgc.googlecode.com/git/img/wsgc_1024_63_noawgn_nofad_trn.4k.png)

_Horizontal lines are artifacts due to the looping of file samples_

# Spectrum in the presence of noise #

## 1023 Hz chip frequency with -27 dB S/N AWGN ##

Now we add AWGN with a S/N of -27 dB in the simulation. Decoding the training sequence is slightly more robust than decoding the message symbols for which -27 dB S/N would give a lot of errors. The following command is used:

`wsgc_generator -s 4096 -c 1023 -C 1020 -t 2.12 -N 4 -I 3 -R 40 -P 2 -L "rcos:13,1.0,2048.0" -n -27 -o test_BPSK.raw --simulate-trn`

**The waterfall looks like:**

![https://wsgc.googlecode.com/git/img/wsgc_4096_1023_-27dB_nofad_trn.wf.png](https://wsgc.googlecode.com/git/img/wsgc_4096_1023_-27dB_nofad_trn.wf.png)

Nothing is visible even with a FFT size of 4096 which corresponds to a bin size of 1 Hz since f<sub>s</sub> is 4096 Hz.

Even at a finer resolution with 16k FFT that is 0.25Hz bins nothing is visible:

![https://wsgc.googlecode.com/git/img/wsgc_4096_1023_-27dB_nofad_trn.16k.png](https://wsgc.googlecode.com/git/img/wsgc_4096_1023_-27dB_nofad_trn.16k.png)

It is not possible to track WSGC signals with the classical FFT and waterfall approach with any practical FFT length. To display something useful one would have to use the results of the pilot correlations in a 3D diagram with time, sequence delay and frequency shift dimensions thus despreading the signal. This would by the way provide an interesting channel sounding profile.

## Decoding the signal ##

At such S/N levels the training sequence can be decoded though as shown in the output of the corresponding wsgc\_test command:

```
wsgc_test -s 4096 -c 1023 -C 1020 -t 2.12 -r 0.0 -N 4 -I 3 -R 40 -P 2 -B 4 -z 2 -F 63 -n -27 -L "rcos:13,1.0,2048.0" --simulate-trn  --cuda

Using CUDA implementation
Using options:
------------- 

Sampling frequency ........:   4096.0
Chip frequency ............:   1023.0
Code length ...............:   1023
Code period ...............:      1.00
Samples/code = FFT size ...:   4096
Code shift ................:   1020
Nb of generated symbols ...:     40
PRNs per pilot averaging ..:      4
Pilot averaging period ....:      4.00
Start PRN number ..........:      3
Tx frequency ..............:      2.12
Initial Rx frequency ......:      0.00
SNR(dB) ...................:    -27.0
Fading Model ..............: None
Modulation ................: BPSK
Nb message symbols ........:     64
Nb service symbols ........:      3
Noise PRN .................:     64 (0)
Pilot PRN gain ............:      0.0 dB
Nb frequency steps ........:     63
Nb frequency sub-steps ....:      8
Frequency range ...........: [-31.0:31.0]
Minor freq step size ......:     0.125
Batch size ................:      4 PRNs
Analysis window size ......:      8 PRNs
Analysis window time ......:      4.00
Pilot PRN 1................:     65 (1)
Pilot PRN 2................:     66 (2)
Message time ..............:    160.00

Lowpass FIR filter model:
Type ......................: Raised cosine
Sampling frequency ........:   4096.0
Cutoff frequency ..........:   2048.0
Nb taps ...................:     13
Rolloff factor ............:      1.00
Tap coefficients:
 0 ( -6): 0.001344
 1 ( -5): -0.000000
 2 ( -4): -0.030310
 3 ( -3): 0.000000
 4 ( -2): 0.280490
 5 ( -1): 0.753080
 6 (  0): 1.000000
 7 (  1): 0.753080
 8 (  2): 0.280490
 9 (  3): 0.000000
10 (  4): -0.030310
11 (  5): -0.000000
12 (  6): 0.001344

Processing options:
 - Using CUDA implementation
 - Simulate synchronization training sequence

Generator polynomials:
G1 = X^10 + X^8 + X^5 + X^1 + 1
G2 = X^10 + X^9 + X^7 + X^6 + X^4 + X^1 + 1

Processing 40 symbols: [   3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   35,   36,   37,   38,   39,   40,   41,   42]

...

Produce signal samples...
Apply lowpass FIR filter
Apply AWGN
Signal power: 3.07801 Noise amplitude: 39.2767
Noise power: 1543.91 S/N (dB): -27.0035
New signal power: 1.00266

...

Training sequence correlation time:  0.324075008 s

Pilot correlation analysis results:
#0: Pilot 1: [([1027:1028] :1027(2) : 4), ([1700] :1700(1) : 1), ([1915] :1915(1) : 1), ([3141] :3141(1) : 1), ([3272] :3272(1) : 1)]
#1: Pilot 1: [([1027:1028] :1027(2) : 3), ([757] :757(1) : 1), ([826] :826(1) : 1), ([969] :969(1) : 1), ([1108] :1108(1) : 1), ([1219] :1219(1) : 1)]
#2: Pilot 1: [([1027:1028] :1028(2) : 3), ([145] :145(2) : 2), ([810] :810(1) : 1), ([1287] :1287(1) : 1), ([2333] :2333(1) : 1)]
#3: Pilot 1: [([1027] :1027(3) : 3), ([1779] :1779(1) : 1), ([2884] :2884(1) : 1), ([3383] :3383(1) : 1), ([3572] :3572(1) : 1), ([3958] :3958(1) : 1)]
#4: Pilot 1: [([1026:1027] :1027(6) : 7), ([1485] :1485(1) : 1)]

...

--- pilot correlation records:
Ai PNi PN# Mag.Max Ph@Mx Ti.. Fi.. delta.F S
 0   0 066 00248.0 -0.54 3141  138  -13.75 N
 0   1 066 00276.4 -0.72 3272  290    5.25 N
 0   2 066 00315.4  1.26 1027  231   -2.12 Y
 0   3 066 00328.6  1.20 1027  231   -2.12 Y
 1   4 066 00298.3  0.84 1028  231   -2.12 Y
 1   5 066 00260.2  1.04 1028  231   -2.12 Y
 1   6 066 00271.4 -2.23 1915  183   -8.12 N
 1   7 066 00225.1 -2.96 1700  136  -14.00 N
 2   8 066 00285.0  1.67  757  186   -7.75 N
 2   9 066 00250.5  0.36 1219  207   -5.12 N
 2  10 066 00241.2 -2.43  969  354   13.25 N
 2  11 066 00228.8 -0.74 1108   39  -26.12 N
 3  12 066 00228.3 -2.47  826  214   -4.25 N
 3  13 066 00260.7  1.07 1027  231   -2.12 Y
 3  14 066 00329.3  1.02 1027  231   -2.12 Y
 3  15 066 00272.5  0.77 1028  231   -2.12 Y
 4  16 066 00363.3  0.95 1028  231   -2.12 N
 4  17 066 00268.0  0.90 1028  231   -2.12 N
 4  18 066 00245.7 -0.39  810  271    2.88 N
 4  19 066 00229.6  0.71 1027  231   -2.12 N
 5  20 066 00224.1  0.96 2333  423   21.88 N
 5  21 066 00240.1 -2.00  145  172   -9.50 N
 5  22 066 00264.7 -2.15  145  172   -9.50 N
 5  23 066 00254.5 -2.81 1287  117  -16.38 N
 6  24 066 00239.6 -2.09 3572  497   31.12 N
 6  25 066 00219.8 -2.24 1779  110  -17.25 N
 6  26 066 00250.1  1.94 3383  349   12.62 N
 6  27 066 00227.8 -1.17 2884  278    3.75 N
 7  28 066 00297.6  0.40 1027  231   -2.12 Y
 7  29 066 00284.1 -1.48 3958  117  -16.38 N
 7  30 066 00276.5  0.20 1027  231   -2.12 Y
 7  31 066 00247.3  0.14 1027  231   -2.12 Y
 8  32 066 00237.6 -2.25 1485  214   -4.25 N
 8  33 066 00288.5 -0.30 1027  231   -2.12 Y
 8  34 066 00283.3 -0.31 1027  231   -2.12 Y
 8  35 066 00450.8 -0.06 1027  231   -2.12 Y
 9  36 066 00347.5 -0.01 1027  231   -2.12 Y
 9  37 066 00431.7  0.09 1027  231   -2.12 Y
 9  38 066 00375.9  0.25 1027  231   -2.12 Y
 9  39 066 00322.9  0.20 1026  231   -2.12 Y
--- training correlation records:
PNi Ai PN# Mag.Max Max/Avg MS Ti..  S
  0  0 000 00000.0 - N/A -  N 3141  N
001 01 000 00000.0 - N/A -  N 3272  N
002 02 000 00000.0 - N/A -  N 1027  Y
003 03 000 00000.0 - N/A -  N 1027  Y
004 04 000 00000.0 - N/A -  N 1028  Y
005 05 000 00000.0 - N/A -  N 1028  Y
006 06 000 00000.0 - N/A -  N 1915  N
007 07 003 00157.1 004.908  Y 1700  N
008 00 000 00000.0 - N/A -  N 0757  N
009 01 000 00000.0 - N/A -  N 1219  N
010 02 000 00000.0 - N/A -  N 0969  N
011 03 000 00000.0 - N/A -  N 1108  N
012 04 000 00000.0 - N/A -  N 0826  N
013 05 000 00000.0 - N/A -  N 1027  Y
014 06 000 00000.0 - N/A -  N 1027  Y
015 07 011 00135.7 006.068  Y 1028  Y
024 00 000 00000.0 - N/A -  N 3572  N
025 01 000 00000.0 - N/A -  N 1779  N
026 02 000 00000.0 - N/A -  N 3383  N
027 03 000 00000.0 - N/A -  N 2884  N
028 04 000 00000.0 - N/A -  N 1027  Y
029 05 000 00000.0 - N/A -  N 3958  N
030 06 000 00000.0 - N/A -  N 1027  Y
031 07 027 00143.9 007.195  Y 1027  Y
032 00 000 00000.0 - N/A -  N 1485  N
033 01 000 00000.0 - N/A -  N 1027  Y
034 02 000 00000.0 - N/A -  N 1027  Y
035 03 000 00000.0 - N/A -  N 1027  Y
036 04 000 00000.0 - N/A -  N 1027  Y
037 05 000 00000.0 - N/A -  N 1027  Y
038 06 000 00000.0 - N/A -  N 1027  Y
039 07 032 00103.0 002.289  Y 1026  Y

```

_The CUDA implementation samples only at the end of analysis period thus the results are relevant for the last PRN in the analysis window only. That is PRNs with analysis index (Ai) of 7._
_Analysis window size is 8 PRNs_

_Maximum of correlation values is scaled to give a reasonable display. The maximum to average is not_

_Although the sequence is shifted by 1020 samples at the start the correlation peak is found at a value around 1027 samples. This is an effect of the additional FIR filter, 7 is about half the number of taps (13)._

_A similar command using the host implementation runs in 3.51377741 s. Compared to the CUDA execution of 0.324075008 s this is a factor of about 10. This factor is observed also for the message correlation_

The whole sequence started at PRN #3 and this is the number found at end of first analysis window. The second analysis window starts at PRN #11 (3+8 = 11). The next at #19 (11+8 = 19), this window is zapped because of poor pilot correlation and no corresponding records are created. Then  at #27 (19+8 = 27). Then at #35 (27+8 = 35), this one is wrongly decoded (32 instead of 35) but in this case the maximum to average (Max/Avg) of correlations is noticeably lower and result can be ruled out. This number is a measure of how well the correlation peak stands out.

If enough samples can be collected the results can be checked against each other as the PRN number found should normally progress by the analysis window size (here 8, starting at 3 this yields the following sequence: 3,11,19,27...)

## 63 Hz chip frequency with -20 dB S/N AWGN ##

One may think the spectrum would become visible in the noise due to its narrower bandwidth. Alas, with a S/N of -20 dB in 1024 Hz bandwidth which correspond to even slightly better conditions than with the previous test it is still invisible as shown by this waterfall with 0.25 Hz bins. Generator command is:

`wsgc_generator -s 1024 -c 63 -m 6 -M 48 -C 50 -t 2.12 -N 4 -I 1 -R 36 -P 2 -L "rcos:21,1.0,128.0" -n -20 -o test_BPSK.raw --simulate-trn`

![https://wsgc.googlecode.com/git/img/wsgc_1024_63_-20dB_nofad_trn.4k.png](https://wsgc.googlecode.com/git/img/wsgc_1024_63_-20dB_nofad_trn.4k.png)

_One can see some patterns in the noise however bear in mind that the file has a limited length and is wrapped around. Also noise is generated by software and is actually pseudo-noise that may not have real noise characteristics_

Test program command is:

`wsgc_test -s 1024 -c 63 -m 6 -C 50 -t 2.12 -r 0.0 -N 4 -I 3 -R 36 -P 2 -B 4 -z 2 -F 63 -M 48 -n -20 --simulate-trn --cuda`