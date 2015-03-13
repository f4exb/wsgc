

# Introduction #

In this page simulation tests using only AWGN are presented.


# Sequence length #

The influence of sequence length is examined here. Results are listed by genrator's maximum power of two m. The sequence length is 2<sup>m</sup>-1.

We try to find the limit SNR for nearly 100% decodes. There is no strict protocol at the moment. Nearly 100% means 1 or 2 faults out of 16 symbols not "very" often that is about 15% of the time.

In most cases we assume a PRN period of 1s. That is f<sub>c</sub> is equal to sequence length in Hz.

PRN integration factor is 4 that is there are 4 PRNs per symbol. A sequence of 16 symbols is produced.

## m=6 ##

This one was a bit difficult to establish it seems not very stable. However for the f<sub>s</sub> listed it was OK-ish.

Typical command:<br>
<code>wsgc_test -s 2048 -c 63 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 63 -n -20 -m 6 -M 58 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 2048 </td><td> 63 </td><td> 16 </td><td> 12 </td><td> -21 </td><td> -9 </td></tr>
<tr><td> 1024 </td><td> 63 </td><td> 8 </td><td> 9 </td><td> -18 </td><td> -9 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth is about <b>-9 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 63 </td><td> 21 </td><td> -153 </td><td> -162 </td></tr></tbody></table>

<h2>m=7</h2>

Here the lowpass filter was active. Its presence has no noticeable influence on the results.<br>
<br>
Typical command:<br>
<code>/wsgc_test -s 1024 -c 127 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 63 -n -18 -L rcos:13,1.0,256.0 -m 7 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 512 </td><td> 127 </td><td> 4 </td><td> 6 </td><td> -17 </td><td> -11 </td></tr>
<tr><td> 1024 </td><td> 127 </td><td> 8 </td><td> 9 </td><td> -20 </td><td> -11 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth: <b>-11 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 127 </td><td> 24 </td><td> -150 </td><td> -161 </td></tr></tbody></table>

<h2>m=8</h2>

Typical command:<br>
<code>/wsgc_test -s 1024 -c 255 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 63 -n -17 -m 8 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 1024 </td><td> 255 </td><td> 4 </td><td> 6 </td><td> -20 </td><td> -14 </td></tr>
<tr><td> 2048 </td><td> 255 </td><td> 8 </td><td> 9 </td><td> -23 </td><td> -14 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth: <b>-14 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 255 </td><td> 27 </td><td> -147 </td><td> -161 </td></tr></tbody></table>

<h2>m=9</h2>

Typical command:<br>
<code>wsgc_test -s 4096 -c 511 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 31 -n -24 -m 9 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 2048 </td><td> 511 </td><td> 2 </td><td> 3 </td><td> -21 </td><td> -18 </td></tr>
<tr><td> 4096 </td><td> 511 </td><td> 4 </td><td> 6 </td><td> -24 </td><td> -18 </td></tr>
<tr><td> 8192 </td><td> 511 </td><td> 8 </td><td> 9 </td><td> -27 </td><td> -18 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth: <b>-18 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 511 </td><td> 30 </td><td> -144 </td><td> -162 </td></tr></tbody></table>

<h2>m=10</h2>

Typical command:<br>
<code>wsgc_test -s 4096 -c 1023 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 63 -n -24 -L rcos:13,1.0,2048.0 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 4096 </td><td> 1023 </td><td> 2 </td><td> 3 </td><td> -24 </td><td> -21 </td></tr>
<tr><td> 8192 </td><td> 1023 </td><td> 4 </td><td> 6 </td><td> -27 </td><td> -21 </td></tr>
<tr><td> 16384 </td><td> 1023 </td><td> 8 </td><td> 9 </td><td> -30 </td><td> -21 </td></tr>
<tr><td> 8192 </td><td> 2046 </td><td> 2 </td><td> 3 </td><td> -24 </td><td> -21 </td></tr>
<tr><td> 4096 </td><td> 511.5 </td><td> 4 </td><td> 6 </td><td> -27 </td><td> -21 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth: <b>-21 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 511.5 </td><td> 30 </td><td> -144 </td><td> -165 </td></tr>
<tr><td> 1023 </td><td> 33 </td><td> -141 </td><td> -162 </td></tr>
<tr><td> 2046 </td><td> 36 </td><td> -138 </td><td> -159 </td></tr></tbody></table>

<h2>m=11</h2>

Typical command:<br>
<code>wsgc_test -s 8192 -c 2047 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 31 -n -24 -m 11 --cuda --simulate-sync</code>

<table><thead><th> <b>f<sub>s</sub></b> </th><th> <b>f<sub>c</sub></b> </th><th> <b>f<sub>s</sub>/2f<sub>c</sub></b> </th><th> <b>adjust</b><br>SNR(dB)</th><th> <b>SNR<sub>lim</sub>(dB)</b> </th><th> <b>SNR<sub>lim</sub>(dB)</b><br>in 2f<sub>c</sub></th></thead><tbody>
<tr><td> 8192 </td><td> 2047 </td><td> 2 </td><td> 3 </td><td> -27 </td><td> -24 </td></tr>
<tr><td> 16384 </td><td> 2047 </td><td> 4 </td><td> 6 </td><td> -30 </td><td> -24 </td></tr></tbody></table>

SNR<sub>lim</sub> in bandwidth: <b>-24 dB</b>.<br>
<br>
Link budget:<br>
<br>
<table><thead><th> <b>f<sub>c</sub>(Hz)</b> </th><th> <b>2f<sub>c</sub>(dB<sub>Hz</sub>)</b> </th><th> <b>N(dBm)</b> </th><th> <b>S(dBm)</b> </th></thead><tbody>
<tr><td> 2047 </td><td> 36 </td><td> -138 </td><td> -159 </td></tr></tbody></table>

<h2>Summary sheets</h2>

<h3>Minimum SNR in signal bandwidth vs sequence length (m)</h3>

<img src='https://wsgc.googlecode.com/git/img/SNR_vs_m.png' />

The trend is nearly regular and proportional to sequence length. Only m=7 and m=8 are slightly off the line only by 1dB.<br>
<br>
<h3>Minimum required power vs sequence length (m)</h3>

<img src='https://wsgc.googlecode.com/git/img/MinPwr_vs_m.png' />

Assuming only thermal noise: -174 dBm/Hz.<br>
<br>
For a constant PRN period length (here 1s) the minimum required power does not vary very much. It is at about -161 / -162 dBm. Longer sequences degrade rapidly because they need a longer PRN period to exploit their advantages. Other sequence lengths are at par or nearly at par because the reduction in bandwidth compensates for the 3dB difference in the limit SNR at every smaller m and in all cases 1s is long enough to have f<sub>c</sub> in Hz equal to the length of the PRN.<br>
<br>
Comparatively JT65B does slightly better at -163 dBm: near 100% decode is obtained at a SNR of -23dB in 2500 Hz according to the chart at p.15 of <a href='http://physics.princeton.edu/pulsar/K1JT/JT65.pdf'>JT65 Communication protocol</a>. 2500Hz is 34 dB<sub>Hz</sub> so the noise power is -174+34 = -140 dBm. Hence the signal power necessary for near 100% decode is -140-23 = -163 dBm. However there is no forward error correction (FEC) with WSGC in this simulation.