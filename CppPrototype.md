

# Introduction #

The prototype written in C++ aims at validating algorithms using locally produced input test samples. Some fading and Gaussian noise can be added to simulate a real radio channel communication.

# Design #
## Overview ##
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_overview.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_overview.png)
## Source generator ##
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_source.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_source.png)
## Fading model ##
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_fading_model.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_fading_model.png)
## Correlator for frequency sensitive modulations (BPSK) ##
![https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_correlator.png](https://wsgc.googlecode.com/git/img/wsgc_cpp_proto_correlator.png)
## Decision box ##

# Program documentation #

Please have a look at the [Doxygen documentation](http://f4exb.free.fr/wsgc/proto_cpp/doc/html/index.html)

# Program usage #

The following binaries are available:
  * wsgc\_test : This is the test program properly speaking. It allows exercising the training sequence and message sequence decoding.
  * wsgc\_generator: This is only the signal generation part and produces an output file of interleaved I and Q floating point samples
  * wsgc\_cuda: This was used to test the CUDA algorithms

## Decoding test: wsgc\_test ##

### wsgc\_test parameters ###

| **Short option** | **Long option** | **Data type** | **Comments** | **Default** |
|:-----------------|:----------------|:--------------|:-------------|:------------|
|  | --help | flag | Displays on line help | false |
|  | --cuda | flag | Triggers CUDA implementation of some (critical) algorithms.<br>Effective only if the program is compiled with CUDA option <table><thead><th> false </th></thead><tbody>
<tr><td>  </td><td> --simulate-trn </td><td> flag </td><td> Triggers the simulation of the training sequence processing. Does message processing if false. </td><td> false </td></tr>
<tr><td> -s </td><td> --f-sampling </td><td> float </td><td> Main sampling frequency </td><td> 4096.0 </td></tr>
<tr><td> -c </td><td> --f-chip </td><td> float </td><td> Chip frequency </td><td> 1023.0 </td></tr>
<tr><td> -C </td><td> --code-shift </td><td> int </td><td> PRN code shift at origin in number of samples </td><td> 1020 </td></tr>
<tr><td> -N </td><td> --nb-prn-per-symbol </td><td> int </td><td> Number of PRN code sequences per symbol </td><td> 4 </td></tr>
<tr><td> -I </td><td> --prn-shift </td><td> int </td><td><i>message simulation</i>: Start simulation at PRN shift in symbol<br><i>training simulation</i>: Start simulation at PRN position in sequence </td><td> 0 </td></tr>
<tr><td> -t </td><td> --f-tx </td><td> float </td><td> Transmission frequency (within baseband) </td><td> 100.0 </td></tr>
<tr><td> -r </td><td> --f-rx </td><td> float </td><td> Initial reception frequency (within baseband) </td><td> 100.1 </td></tr>
<tr><td> -n </td><td> --snr </td><td> float </td><td> Signal to noise ratio in dB </td><td> -24.0 </td></tr>
<tr><td> -m </td><td> --nb-stages </td><td> int </td><td> Number of stages of the PRN generator LFSRs. This gives 2N-1 chips </td><td> 10 </td></tr>
<tr><td> -g </td><td> --g1-poly </td><td> string </td><td> LFSR generator 1 comma separated list of polynomial powers except N and 0 </td><td> "8,5,1" </td></tr>
<tr><td> -G </td><td> --g2-poly </td><td> string </td><td> LFSR generator 2 comma separated list of polynomial powers except N and 0 </td><td> "9,7,6,4,1" </td></tr>
<tr><td> -M </td><td> --nb-message-symbols </td><td> int </td><td> Number of message symbols </td><td> 64 </td></tr>
<tr><td> -S </td><td> --nb-service-symbols </td><td> int </td><td> Number of service symbols. Includes noise and pilot PRNs. </td><td> 3 </td></tr>
<tr><td> -Y </td><td> --nb-training-symbols </td><td> int </td><td> Number of training symbols, used with --simulate-trn and unpiloted schemes </td><td> 0 </td></tr>
<tr><td> -P </td><td> --pilot-prns </td><td> int </td><td> Number of pilot PRNs.<br>Use 1 for piloted message correlation.<br>Use 2 for training sequence correlation.</td><td> 1 </td></tr>
<tr><td> -A </td><td> --pilot-gain-db </td><td> float </td><td> Gain of the pilot PRN(s) vs message PRNs in dB. Usually zero or negative. </td><td> 0.0 </td></tr>
<tr><td> -F </td><td> --df-steps </td><td> int </td><td> Number of frequency steps to explore.<br>Step sized by FFT size and sampling frequency.<br>Should be an odd number >= 3 so as to have center frequency centered. </td><td> 31 </td></tr>
<tr><td> -U </td><td> --df-sub-steps </td><td> int </td><td> Number of frequency sub-steps of the above frequency step to explore </td><td> 8 </td></tr>
<tr><td> -B </td><td> --batch-size </td><td> int </td><td> Batch size in number of PRNs for some algorithms implemented with batch processing.<br>These are the pilot and "piloted" message processing. </td><td> 4 </td></tr>
<tr><td> -p </td><td> --prns </td><td> string </td><td> Given symbol PRN numbers in a comma separated list.<br><i>message simulation:</i> Overrides -R option<br><i>training simulation:</i> Ignored </td><td> (null) </td></tr>
<tr><td> -R </td><td> --random-prns </td><td> int </td><td> <i>message simulation:</i> Number of random symbols to generate.<br><i>training simulation:</i> Number of PRNs to generate in the training sequence </td><td> 4 </td></tr>
<tr><td> -f </td><td> --fading-model </td><td> string </td><td> Fading model parameters (see next paragraph).<br>No fading (flat signal) if not specified </td><td> (null) </td></tr>
<tr><td> -d </td><td> --modulation-scheme </td><td> string </td><td> BPSK: Binary Phase Shift Keying. Bits modulate the signal at phase 0 (1) and 180 degrees (0 or -1). Yields piloted correlation scheme.<br>OOK: On-Off Keying. The carrier is either present for bit 1 or absent for bit 0. Yields unpiloted correlation scheme.<br>DBPSK: Differential BPSK phase of chip is used as reference for next chip. Otherwise BPSK. Yields unpiloted correlation scheme.<br>MFSK: Multiple Frequency Shift Keying (similar to JT65) takes 3 arguments (see next) </td><td> BPSK </td></tr>
<tr><td> -z </td><td> --analysis-window-size </td><td> int </td><td> Pilot analysis window size in number of symbols.<br>A symbol is comprise of the repetition of several PRNs (see -N option). </td><td> 4 </td></tr>
<tr><td> -L </td><td> --fir-filter-model </td><td> string </td><td> Use a lowpass FIR filter with given characteristics for pulse shaping (see next paragraph).<br>No filtering if not specified </td><td> (null) </td></tr>
<tr><td> -y </td><td> --gpu-affinity </td><td> int </td><td> Specify on which GPU to run in CUDA mode </td><td> (default attribution) </td></tr>
<tr><td> -j </td><td> --source-coding </td><td> string </td><td> Applies source coding to a given text message to produce message bytes (see parameters next). Overrides -R, -p, -M and possible Reed-Solomon encoding options. Other parameters like -m are checked for consistency </td><td> (null) </td></tr>
<tr><td> -O </td><td>--fec-option </td><td> string </td><td> Applies Forward Error Correction using external library (see parameters next) </td><td> (null) </td></tr></tbody></table>

<h3>Default LFSR generator polynomials</h3>

By default we use polynomials taken among the so called "preferred m-sequences" that give best ratio between correlation peak and maximum positive and negative peaks for uncorrelated sequences. For each number of stages (or maximum power of two of the polynomial) we have:<br>
<br>
<table><thead><th> <b>Nb stages</b><br>m (2<sup>m</sup>)</th><th> <b>G1</b> </th><th> <b>G2</b> </th><th> <b>Nb PRNs</b><br>2<sup>m</sup>+1</th><th> <b>PRN len</b><br>2<sup>m</sup>-1</th><th> <b>Autocorr</b><br>Hi</th><th> <b>Autocorr</b><br>Low</th><th> <b>Crosscorr</b><br>Hi</th><th> <b>Crosscorr</b><br>Low</th><th> <b>Bound</b><br>(%)</th></thead><tbody>
<tr><td> 5 (32) </td><td> X<sup>5</sup> + X<sup>3</sup> + 1 </td><td> X<sup>5</sup> + X<sup>4</sup> + X<sup>3</sup> + X<sup>2</sup> + 1 </td><td> 33 </td><td> 31 </td><td> 31 </td><td> -1 </td><td> 7 </td><td> -9 </td><td> -29 </td></tr>
<tr><td> 6 (64)</td><td> X<sup>6</sup> + X<sup>1</sup> + 1 </td><td> X<sup>6</sup> + X<sup>5</sup> + X<sup>2</sup> + X<sup>1</sup> + 1 </td><td> 65 </td><td> 63 </td><td> 63 </td><td> -1 </td><td> 15 </td><td> -17 </td><td> -27 </td></tr>
<tr><td> 7 (128)</td><td> X<sup>7</sup> + X<sup>3</sup> + 1 </td><td> X<sup>7</sup> + X<sup>3</sup> + X<sup>2</sup> + X<sup>1</sup> + 1 </td><td> 129 </td><td> 127 </td><td> 127 </td><td> -1 </td><td> 15 </td><td> -17 </td><td> -13 </td></tr>
<tr><td> 8 (256) </td><td> X<sup>8</sup> + X<sup>7</sup> + X<sup>6</sup> + X<sup>1</sup> + 1 </td><td> X<sup>8</sup> + X<sup>7</sup> + X<sup>6</sup> + X<sup>5</sup> + X<sup>4</sup> + X<sup>2</sup> + 1 </td><td> 257 </td><td> 255 </td><td> 255 </td><td> -1 </td><td> 31 </td><td> -17 </td><td> 12 </td></tr>
<tr><td> 9 (512) </td><td> X<sup>9</sup> + X<sup>4</sup> + 1 </td><td> X<sup>9</sup> + X<sup>6</sup> + X<sup>4</sup> + X<sup>3</sup> + 1 </td><td> 513 </td><td> 511 </td><td> 511 </td><td> -1 </td><td> 31 </td><td> -33 </td><td> -6 </td></tr>
<tr><td> 10 (1024) </td><td> X<sup>10</sup> + X<sup>8</sup> + X<sup>5</sup> + X<sup>1</sup> + 1 </td><td> X<sup>10</sup> + X<sup>9</sup> + X<sup>7</sup> + X<sup>6</sup> + X<sup>4</sup> + X<sup>1</sup> + 1 </td><td> 1025 </td><td> 1023 </td><td> 1023 </td><td> -1 </td><td> 63 </td><td> -65 </td><td> -6 </td></tr>
<tr><td> 11 (2048) </td><td> X<sup>11</sup> + X<sup>2</sup> + 1 </td><td> X<sup>11</sup> + X<sup>8</sup> + X<sup>5</sup> + X<sup>2</sup> + 1 </td><td> 2049 </td><td> 2047 </td><td> 2047 </td><td> -1 </td><td> 63 </td><td> -65 </td><td> -3 </td></tr></tbody></table>

<h3>Fading options</h3>

The general format for fading options is:<br>
<code>&lt;fading model character code&gt;:&lt;fading model parameters&gt;</code>
The fading model character code is:<br>
<ul><li>C: Clarke's model<br>
</li><li>W: Watterson's model</li></ul>

<h4>Clarke's model</h4>

Format is <code>C:&lt;nb paths&gt;,&lt;frequency spread (Hz)&gt;</code>

For a more random like fading envelope you will need at least 8 paths.<br>
<br>
<h4>Watterson's model</h4>

Format is <code>W:&lt;path params&gt;/&lt;path params&gt;/...</code>
Path params are: <code>&lt;delay (s)&gt;,&lt;amplitude&gt;,&lt;f Doppler (Hz)&gt;,&lt;f offset (Hz)&gt;</code>

To simulate a ionospheric channel you would set at least a pair of paths with same amplitude factor with one delayed with respect to the other and identical Doppler frequency spreading. The offset frequency can be set to a sub-Hertz value.<br>
<br>
You may simulate a tropospheric channel pretty much the same way but the delay is going to be considerably smaller. You can also omit the frequency offset. You may specify a lower amplitude for the delayed channel(s) to simulate a main path with secondary echoes.<br>
<br>
You can check typical values for delay (time spreading) and Doppler frequency spreading <a href='http://code.google.com/p/wsgc/wiki/FadingModels#Typical_fading_characteristics_for_VHF_and_up'>here</a>

<h3>Lowpass FIR filter for pulse shaping</h3>

A BPSK modulated signal spectrum although infinitely decreasing is infinite. It is often interesting if not required to limit the bandwidth of the BPSK signal to just what is necessary for proper decoding. The essential part of the signal for that matter is in the main lobe around zero frequency. It is practical to choose a FIR filter and the type of filter is carefully chosen to keep good decoding properties (low inter-symbol interference). One of such filters is the raised cosine filter that is the only one implemented here for the moment.<br>
<br>
General format is: <code>&lt;Filter type&gt;:&lt;comma separated parameters&gt;</code>

The number of taps and cutoff frequency are common to all types of filter and therefore always specified in the parameters. However their exact placement and influence may vary.<br>
<br>
The output of the programs show the filter characteristics to help tuning.<br>
<br>
<h4>Raised cosine</h4>

Its essential characteristics are the number of taps and the roll-off factor.<br>
<br>
Format is: <code>rcos:&lt;nb taps&gt;,&lt;roll-off factor&gt;,&lt;cutoff frequency&gt;</code>

<ul><li>Cutoff frequency must not exceed half the sampling frequency<br>
</li><li>For complex signals as this is the case the cutoff frequency must be at least twice the bit rate (that is the chip rate).<br>
</li><li>Roll-off factor: the most aggressive shape is obtained with a roll-off factor of 1.0 while a roll-off of 0.0 will make the transition very smooth. For this application a roll-off of 1.0 is perfectly convenient.<br>
</li><li>The number of taps relate to the number of chips that will be involved. You can choose it so that edge values are close to zero to allow a soft transition.</li></ul>


<h3>Modulation options</h3>

<h4>BPSK, OOK, DBPSK</h4>

Take no specific options<br>
<br>
<h4>MFSK</h4>

Takes three options after a colon separator i.e. <code> -d MFSK:&lt;log2 bandwidth&gt;,&lt;log2 symbol time&gt;,&lt;zero slot shift&gt; </code>
<ul><li><b>log2 bandwidth</b>: log2 of frequency tone separation hence symbol bandwidth in Hz<br>
</li><li><b>log2 symbol time</b>: log2 of one symbol (or tone) time in seconds<br>
</li><li><b>zero slot shift</b>: shift of symbol 0 in number of frequency slots with respect to slot at center frequency.</li></ul>

The system will adjust the FFT size so that one FFT bin width equals the symbol bandwidth.<br>
<br>
Ex: <code>-s 4096 -d MFSK:3,0,-32 -M 64 -R 63</code> yields:<br>
<pre><code>Modulation ................: MFSK<br>
Sampling frequency ........:   4096.0<br>
Symbol bandwidth ..........:      8.0<br>
Zero frequency ............:   -256.0<br>
Symbol time ...............:      1.000<br>
FFT size ..................:    512<br>
Zero FFT slot .............:    -32<br>
Nb FFT per symbol .........:      8<br>
Nb message symbols ........:     64<br>
Nb service symbols ........:      3<br>
Nb training symbols .......:      0<br>
Max symbol frequency.......:    280.00<br>
Nb of generated symbols ...:     63<br>
Message time ..............:     63.00<br>
Tx shift frequency ........:      0.00<br>
</code></pre>

Note on MFSK: the receiving frequency shift is not implemented, therefore transmissions must be centered on 0 Hz i.e. use option <code>-t 0.0</code>

<h3>Source coding options</h3>

<code>-j &lt;source codec code&gt;:&lt;source codec parameters&gt;:&lt;textual message&gt;</code></br></br>
The source codec code is:<br>
<ul><li><b>JT65</b>: This is the original JT65 encoding and takes no parameters. It implies -M 64 and RS(63,12) code if Reed-Solomon is used in which case it generates 63 symbols. The number of GC generator stages must be higher than 6 (-m parameter).<br>
</li><li><b>JT257</b>: This is a variant of the original JT65 packing the 72 bits into 9 symbols of 8 bits. It implies -M 256 and RS(255,9) if Reed-Solomon is used in which case it generates 255 symbols.. The number of GC generator stages must be higher than 8 (-m parameter).<br>
</li><li><b>JTCC</b>: The 72 bits are not packed that is they are organized individually as 72 1 bit symbols. This is normally used with rate 1/n convolutional coding such as the rate 1/2 k=32 Layland-Lushbaugh code of JT4.</li></ul>

<h3>FEC options</h3>
<h4>Reed-Solomon using RSSoft library</h4>
<code>-O RS/q,k,M,r,i:&lt;increment mode&gt;,&lt;decoding mode&gt;(,&lt;optional argument&gt;)</code></br></br>Options are:<br>
<ul><li><b>q</b>: n = 2^q-1 for RS(n,k)<br>
</li><li><b>k</b>: k for RS(n,k)<br>
</li><li><b>M</b>: Initial global multiplicity for multiplicity matrix<br>
</li><li><b>r</b>: If no solution is found the soft-decision algorithm is relaunched with a larger global multiplicity of the multiplicity matrix. This is the maximum number of retries<br>
</li><li><b>i</b>: Global multiplicity base increment<br>
</li><li>global multiplicity <b>increment mode</b> is one of the following:<br>
<ul><li><b>arith</b>: next multiplicity is multiplicity plus base increment<br>
</li><li><b>iarith</b>: next multiplicity is multiplicity plus increment that is previous increment plus base increment<br>
</li><li><b>geom</b>: next multiplicity is multiplicity times base increment<br>
</li><li><b>igeom</b>: next multiplicity is multiplicity times increment that is previous increment times base increment<br>
</li></ul></li><li><b>decoding mode</b> is one of the following:<br>
<ul><li><b>all</b>: returns all messages with their probability score running for the maximum number of retries (full scan)<br>
</li><li><b>full</b>: returns unique messages with their best probability score running for the maximum number of retries (full scan)<br>
</li><li><b>best</b>: returns only the message with best probability score (this implies  it runs for the maximum of retries)<br>
</li><li><b>first</b>: returns only the first message found<br>
</li><li><b>regex</b>: runs for the maximum number of retries and until it finds the message matching the given regular expression. In this case the regular expression is given next after the comma separator.<br>
</li><li><b>match</b>: runs for the maximum number of retries and until it finds the message matching exactly the given string. In this case the matching string is given next after the comma separator.<br>
</li><li><b>binmatch</b>: runs for the maximum number of retries and until it finds the message matching exactly the message that was sent. It is meant to be used in a pure testing context.<br>
</li><li><b>relthr</b>: runs for the maximum number of retries and until it finds the message having a probability score (or reliability) greater than the specified threshold. In this case the reliability threshold as a decimal value is given next after the comma separator.<br>
<h4>Convolutional coding using CCSoft library</h4>
<code>-O CC/&lt;constraint lengths&gt;/&lt;generator polynomials&gt;/&lt;decoding parameters&gt;</code></br></br>Options are:<br>
</li></ul></li><li><b>constraint lengths:</b> comma separated list of constraint (actually register) lengths. there is one constraint per input symbol bit (k).<br>
</li><li><b>generator polynomials:</b> colon separated list of lists of binary representation of generator polynomials.<br>
<ul><li>there is one list per input symbol bit (k)<br>
</li><li>each list has one polynomial per output symbol bit (n)<br>
</li></ul></li><li><b>decoding parameters:</b> <code>&lt;</code>algorithm code<code>&gt;</code>,<code>&lt;</code>interleaving option<code>&gt;</code>,<code>&lt;</code>decoding option<code>&gt;</code>:<code>&lt;</code>decoding arguments comma separated<code>&gt;</code> . All arguments have defaults and are not mandatory<br>
<ul><li><b>interleaving option</b>:<br>
<ul><li><b>ni</b>: No interleaving<br>
</li><li><b>(empty)</b>: Default i.e. index bit reversal interleaving<br>
</li></ul></li><li><b>decoding option</b>:<br>
<ul><li><b>(empty)</b>: No retry mechanism<br>
</li><li><b>regex</b>: Retry with regular expression match<br>
</li><li><b>matchstr</b>: Retry with exact textual message match<br>
</li><li><b>matchmsg</b>: Retry with exact binary message match (for benckmarks)<br>
</li></ul></li><li><b>algorithm code</b>:<br>
<ul><li><b>stack:</b> Stack algorithm. Arguments are: <code>&lt;</code>bias<code>&gt;</code>,<code>&lt;</code>node limit<code>&gt;</code>,<code>&lt;</code>metric limit<code>&gt;</code>,<code>&lt;</code>maximum nb of retries<code>&gt;</code>,<code>&lt;</code>edge metric bias decrement<code>&gt;</code>
<ul><li><b>bias</b>: Edge metric bias. Floating point default 0.0<br>
</li><li><b>node limit</b>: Maximum number of nodes in the code tree. Unsigned integer default 0 (unused)<br>
</li><li><b>metric limit</b>: Stops when the path metric falls below this value. Floating point no default (unused)<br>
</li><li><b>maximum number of retries</b>: When a retry decoding option is used this is the maximum number of retries<br>
</li><li><b>edge metric bias decrement</b>: At each retry the edge bias is decremented by this value<br>
</li></ul></li><li><b>fano:</b> Fano algorithm. Arguments are: <code>&lt;</code>bias<code>&gt;</code>,<init threshold>,<delta threshold>,<cache size>,<delta init threshold>,<node limit>,<threshold limit>,<code>&lt;</code>maximum nb of retries<code>&gt;</code>,<code>&lt;</code>edge metric bias decrement<code>&gt;</code>
<ul><li><b>bias</b>: Edge metric bias. Floating point default 0.0<br>
</li><li><b>init threshold</b>: Initial path metric threshold. Floating point default 0.0<br>
</li><li><b>delta threshold</b>: Amount by which the path metric threshold is incrementend (tightened) or decremented (loosened). Floating point default 1.0<br>
</li><li><b>cache size</b>: Maximum number of nodes kept in memory at once. Unsigned integer default 0 (cache unused)<br>
</li><li><b>delta init threshold:</b> Amount added to threshold to restart process when a loop condition is detected. To be effective this number has to be negative (decrement threshold). Floating point default 0.0 (inactive)<br>
</li><li><b>node limit</b>: Maximum number of nodes in the code tree. Unsigned integer default 0 (unused)<br>
</li><li><b>threshold limit</b>: Stops when the current path metric threshold falls below this value. Floating point no default (unused)<br>
</li><li><b>maximum number of retries</b>: When a retry decoding option is used this is the maximum number of retries<br>
</li><li><b>edge metric bias decrement</b>: At each retry the edge bias is decremented by this value</li></ul></li></ul></li></ul></li></ul>

<u>Retry mechanism</u>: the decoding options regex, matchstr and matchmsg implement a retry mechanism. Until the match condition is found or the maximum number of retries is reached, the decoding process is restarted with a lower edge metric bias. Edge metric starts at the given edge metric and is lowered by the edge metric bias decrement at each retry.<br>
<h3>Command line examples</h3>
<ul><li>Message correlation OOK:</br><code>wsgc_test -s 4096 -c 1023 -C 10 -t 0.0 -r 0.0 -N 4 -I 0 -p 1,3,5,7,9,11,13,15 -B 2 -z 4 -n -6 -d OOK -H 0.13,1.1</code>
</li><li>Message correlation DBPSK:</br><code>wsgc_test -s 4096 -c 1023 -C 10 -t 0.0 -r 0.0 -N 4 -I 0 -p 1,3,5,7,9,11,13,15 -B 2 -z 4 -n -13 -d DBPSK -H 0.13,1.1</code>
</li><li>Message correlation BPSK with CUDA:</br><code>wsgc_test -s 4096 -c 1023 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 31 -n -24 --cuda</code>
</li><li>Message correlation MFSK:</br><code>wsgc_test -s 4096 -d MFSK:3,-1,-32 -t 0.0 -M 64 -R 63 -n -24 -H 1.75</code>
</li><li>MFSK with JT65 classic source encoding and Reed-Solomon with RSSoft library:</br><code>wsgc_test -s 4096 -d "MFSK:3,0,-32" -t 0.0 -M 64 -R 16 -n -20 -f "W:0.0,1.0,5.0,0.0" -j "JT65::F6HTJ F4EXB JN33" -O RS/63,12,16,5,2:iarith,regex,.*F4EXB.*</code>
</li><li>MFSK with JT65 1 bit encoding (JT4 like) and Convolutional Coding (CCSoft library) with retry option (matchmsg):</br><code>wsgc_test -s 4096 -d "MFSK:1,0,-32" -t 0.0 -M 2 -j "F6HTJ F4EXB JN33" -n -30.9 -O CC/32/4073739089,3831577671/stack,,matchmsg:-1.5,2000000,-24,5,0.05</code></li></ul>

<h2>Generating samples file: wsgc_generator</h2>

<h3>wsgc_generator parameters</h3>

Parameters are common with wsgc_test however some are irrelevant while only -o (output file) is specific to wsgc_generator. Relevant commands are listed here:<br>
<br>
<table><thead><th> <b>Short option</b> </th><th> <b>Long option</b> </th><th> <b>Data type</b> </th><th> <b>Comments</b> </th><th> <b>Default</b> </th></thead><tbody>
<tr><td>  </td><td> --help </td><td> flag </td><td> Displays on line help </td><td> false </td></tr>
<tr><td>  </td><td> --simulate-trn </td><td> flag </td><td> Triggers the simulation of the training sequence processing </td><td> false </td></tr>
<tr><td> -s </td><td> --f-sampling </td><td> float </td><td> Main sampling frequency </td><td> 4096.0 </td></tr>
<tr><td> -c </td><td> --f-chip </td><td> float </td><td> Chip frequency </td><td> 1023.0 </td></tr>
<tr><td> -C </td><td> --code-shift </td><td> int </td><td> PRN code shift at origin in number of samples </td><td> 1020 </td></tr>
<tr><td> -N </td><td> --nb-prn-per-symbol </td><td> int </td><td> Number of PRN code sequences per symbol </td><td> 4 </td></tr>
<tr><td> -I </td><td> --prn-shift </td><td> int </td><td><i>message simulation</i>: Start simulation at PRN shift in symbol<br><i>training simulation</i>: Start simulation at PRN position in sequence </td><td> 0 </td></tr>
<tr><td> -t </td><td> --f-tx </td><td> float </td><td> Transmission frequency (within baseband) </td><td> 100.0 </td></tr>
<tr><td> -n </td><td> --snr </td><td> float </td><td> Signal to noise ratio in dB </td><td> -24.0 </td></tr>
<tr><td> -m </td><td> --nb-stages </td><td> int </td><td> Number of stages of the PRN generator LFSRs. This gives 2N-1 chips </td><td> 10 </td></tr>
<tr><td> -g </td><td> --g1-poly </td><td> string </td><td> LFSR generator 1 comma separated list of polynomial powers except N and 0 </td><td> "8,5,1" </td></tr>
<tr><td> -G </td><td> --g2-poly </td><td> string </td><td> LFSR generator 2 comma separated list of polynomial powers except N and 0 </td><td> "9,7,6,4,1" </td></tr>
<tr><td> -M </td><td> --nb-message-symbols </td><td> int </td><td> Number of message symbols </td><td> 64 </td></tr>
<tr><td> -S </td><td> --nb-service-symbols </td><td> int </td><td> Number of service symbols. Includes noise and pilot PRNs. </td><td> 3 </td></tr>
<tr><td> -P </td><td> --pilot-prns </td><td> int </td><td> Number of pilot PRNs.<br>Use 1 for piloted message correlation.<br>Use 2 for training sequence correlation.</td><td> 1 </td></tr>
<tr><td> -A </td><td> --pilot-gain-db </td><td> float </td><td> Gain of the pilot PRN(s) vs message PRNs in dB </td><td> 0.0 </td></tr>
<tr><td> -p </td><td> --prns </td><td> string </td><td> Given symbol PRN numbers in a comma separated list.<br><i>message simulation:</i> Overrides -R option<br><i>training simulation:</i> Ignored </td><td> (null) </td></tr>
<tr><td> -R </td><td> --random-prns </td><td> int </td><td> <i>message simulation:</i> Number of random symbols to generate.<br><i>training simulation:</i> Number of PRNs to generate in the training sequence </td><td> 4 </td></tr>
<tr><td> -f </td><td> --fading-model </td><td> string </td><td> Fading model parameters.<br>No fading (flat signal) if not specified </td><td> (null) </td></tr>
<tr><td> -d </td><td> --modulation-scheme </td><td> string </td><td> BPSK: Binary Phase Shift Keying. Bits modulate the signal at phase 0 (1) and 180 degrees (0 or -1).<br>OOK: On-Off Keying. The carrier is either present for bit 1 or absent for bit 0. </td><td> BPSK </td></tr>
<tr><td> -L </td><td> --fir-filter-model </td><td> string </td><td> Use a lowpass FIR filter with given characteristics for pulse shaping.<br>No filtering if not specified </td><td> (null) </td></tr>
<tr><td> -o </td><td> --samples-output-file </td><td> string </td><td> Path to output file. Mandatory parameter. </td><td> (null) </td></tr>
<tr><td> -j </td><td> --source-coding </td><td> string </td><td> Applies source coding to a given text message to produce message bytes (see parameters next). Overrides -R, -p, -M and possible Reed-Solomon encoding options. Other parameters like -m are checked for consistency </td><td> (null) </td></tr>
<tr><td> -O </td><td> --fec-option </td><td> string </td><td> Applies Forward Error Correction using external library. Only encoding parameters are required </td><td> (null) </td></tr>