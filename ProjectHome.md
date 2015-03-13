![http://nutrition-edu.com/ee/terms/images/Gold%20code.gif](http://nutrition-edu.com/ee/terms/images/Gold%20code.gif)



# Bad news #

A recent benchmark (August 2013) proves that plain non-coherent MFSK which a well known implementation is JT65 does much better than the Gold Code sequences correlation scheme of WSGC in all practical situations. So it is clearly not worth the effort running this complex and CPU burner engine to try digging out signals from the noise. WSGC is just an(other) ill-fated attempt at a weak-signal transmission mode.

More information on:
  * [WSGC vs MFSK](https://code.google.com/p/wsgc/wiki/WSGC_vs_MFSK)
  * [Doppler spread with GC correlation](https://code.google.com/p/wsgc/wiki/MFSK_Doppler_spread#GC_correlation)

It is perhaps useful for something may be for channel sounding. There is no point however in developing it any further as a weak signal transmission mode.

However the code has interesting constructs for Cuda and a fading engine to make channel simulations. It also has a practical implementation of the soft-decision Reed Solomon decoder library [RSSoft](https://code.google.com/p/rssoft).

It could evolve as a weak signal testing platform named **Weak Signal Great Circus** in order to keep the acronym. Below is the future logo. Look at this clown with big ears to hear weak signals. Isn't it cute? Looks like me:

![http://clipartist.info/www/colouringbook/B/big_earred_clown_black_white_line_art_coloring_book_colouring-69px.png](http://clipartist.info/www/colouringbook/B/big_earred_clown_black_white_line_art_coloring_book_colouring-69px.png)

# Introduction #

WSGC stands for "Weak Signal communications using Gold Codes". It is primarily
intended for Amateur Radio (Ham Radio) use to experiment with radio
transmission channels with very marginal signal over noise ratio and with
possibly large fading effects both in the time and frequency domain.

Gold Codes are at the heart of the system to convey information. The Gold
Code sequence itself is the symbol transmitted. The technique was inspired
by the GPS system where the spread weak signal is correlated with a
collection of pre-definite Gold Code sequences to produce an identifiable
correlation peak. In the GPS system a Gold Code sequence identifies the
Satellite Vehicle (SV) number. Here this information itself is used as a
symbol in the alphabet of transmitted symbols (in comparison the fleet of
SVs available).

This project was started in September 2011 as an individual effort to
pursue the opportunity to construct a transmission mode inspired by the
GPS transmission system capable of challenging the very extreme limits of weak
signal transmission possibilities. As a hobbyist project its progress is
subject to individual priorities and availability. Until now the idea has
been exposed and discussed within a small circle of fellow Ham radio
operators in the local Radio Club and a software prototype has been
written.

![http://thrust.googlecode.com/hg/doc/thrust_logo.png](http://thrust.googlecode.com/hg/doc/thrust_logo.png)

The software prototype is based on signal simulation to exercise the
different algorithms and validate the proof of concept. A first version in
Python language then another version closer to the final goals was
developed in C++. More recently the most processing demanding parts have
been implemented using GPU acceleration with nVidia's [CUDA/Thrust](http://thrust.github.com/) framework.

At this point of time (December 2012) a working prototype exists proving
it is possible to correctly receive symbols with a SNR of as low as -27dB subject to AWGN
provided strong error correction is put in place. So it is time to start
thinking of sharing this with the global Ham radio community and possibly
at a larger scale if someone thinks this could be useful for some
different purpose. This software prototype is the very reason why this
project has been opened in Google Code. The whole system apart from the
radio front end is expected to be done in software and hosted here. Documentation about the hardware will also be included.

For the moment a document detailing the process of WSGC signals has been
written in French only. Of course an English version is in the plans. The
software however has been commented in English.

If you are not familiar with Spread Spectrum techniques you may want to check this interesting [Introduction article about Spread Spectrum](http://www.sss-mag.com/pdf/Ss_jme_denayer_intro_print.pdf)

# How does this compare to JT65? #

![http://www.nitehawk.com/rasmit/k2uyh_jt_spec.png](http://www.nitehawk.com/rasmit/k2uyh_jt_spec.png)

_JT65 is the digital mode of choice for weak signal work. You can check all the details on [WSJT Home page](http://physics.princeton.edu/pulsar/K1JT/)._

_WSGC uses a spread spectrum technique and the signal cannot be visualized using a FFT waterfall display like this. See: [WSGC spectrum wiki](https://code.google.com/p/wsgc/wiki/Spectrum)_

At this point this is probably the question you are willing to ask as I did myself. For one aspect it does not really compare because of a complete different approach using a spread spectrum technique instead of narrowband MFSK.

Now that a Reed-Solomon decoder is available using RSSoft more precise comparisons can be made. To put things on par the respective sampling frequencies and additive noise bandwidths have to be considered as we wull see in detail in the next paragraphs.

In a different note WSGC appears to be more demanding in terms of computer resources and for speedy versions in terms of chip rate GPU or even multi-GPU support may be required. This rules out autonomous portable work or lightweight operating conditions in general which JT65 is able to support directly. In some cases reception of WSGC signals will only be supported in remote connection with a "supercomputer" server doing the actual decoding work. Thus the final architecture will have to be able to support remote connection. The radio hardware front end also appears to be considerably more complex than a commercial transceiver and a sound card.

To see a preliminary benchmark using only AWGN you can check a [chart showing the minimum required power vs sequence length](https://code.google.com/p/wsgc/wiki/SimulationAWGN#Minimum_required_power_vs_sequence_length_(m))

## MFSK, JT65 and RSSoft study ##

![https://wsgc.googlecode.com/git/img/RSSoft_logo.png](https://wsgc.googlecode.com/git/img/RSSoft_logo.png)

This part does not concern the correlation engine so it is not part of WSGC properly speaking. However a MFSK non-coherent modulation scheme has been added to the palette of modes supported in the WSGC prototype in order to be able to make some comparisons. This MFSK modulation scheme is quite similar to what happens with JT65 but uses the [RSSoft](https://code.google.com/p/rssoft) soft-decision Reed-Solomon decoder instead of the proprietary KVASD. RSSoft has a 0.5 dB performance hit versus KVASD but is Open Source and can accomodate different rates of Reed-Solomon codes than just the RS(63.12).

Therefore a few simulations were made and reported in the Wiki which serve also to test the RSSoft engine. The following case studies are proposed:
  * [RSSoft vs KVASD](https://code.google.com/p/wsgc/wiki/RSSoft_vs_KV): Compares RSSoft and KVASD performance and shows a 0.5 dB disadvantage for RSSoft. It also sets the base of comparison between WSJT/JT65 figures and WSGC prototype figures.
  * [JT modes galore!](https://code.google.com/p/wsgc/wiki/MFSK_JT_modes_galore): Exploring even more variants of JT modes taking advantage of RSSoft flexibility with respect to RS parameters while also tweaking symbol bandwidth and symbol length.
  * [Doppler spread study](https://code.google.com/p/wsgc/wiki/MFSK_Doppler_spread): How do all these modes behave in fading conditions with Doppler spread. Surprise, surprise!


## MFSK vs WSGC study ##

It seems that plain non-coherent MFSK does much better than the correlation scheme of WSGC at least with just AWGN noise. So it is clearly not worth the effort running this complex and CPU burner engine to try digging out signals from the noise.

More information on:
  * [WSGC vs MFSK](https://code.google.com/p/wsgc/wiki/WSGC_vs_MFSK)
  * [Doppler spread with GC correlation](https://code.google.com/p/wsgc/wiki/MFSK_Doppler_spread#GC_correlation)

# Similar Works #

http://www.merl.com/template/image.php?src=projects/images/iruwbdesign.jpg&width=250

The idea of using the spreading code to carry a few bits of the symbols is not entirely new. The originality of WSGC is to use the spreading code exclusively as I have seen with a few simulations that it is very difficult to recover the main modulation in very weak signal conditions.

There is quite a significant amount of academic research in this area. This technique is usually known as "Code Shift Keying" (CSK) or "M-ary Code Shift Keying" (MCSK). This is indeed quite similar to FSK and MFSK where the difference is that the spreading "code" is used in the place of the signal "frequency" as the symbol information. One of the flagship applications is with Ultra Wide Band transmission techniques in order to carry more bits of information on the same channel. In this case the mainstream of bits is transmitted in the original modulated carrier and the spreading code carries a few additional bits.

You can check a very interesting thesis paper published in 1998 here: [Low Complexity Code Shift Keying Over The Wireless Channel](http://dspace.ucalgary.ca/bitstream/1880/26208/1/38623Chow.pdf) by Mable Man Chee Chow at the Univestity of Calgary.

# Roadmap and News #

![http://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Road_in_Norway-1.jpg/320px-Road_in_Norway-1.jpg](http://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Road_in_Norway-1.jpg/320px-Road_in_Norway-1.jpg)

## Python Prototype ##

Done but deprecated hence it is not stored here

## C++ Prototype ##
[C++ Prototype wiki](http://code.google.com/p/wsgc/wiki/CppPrototype)
[Doxygen documentation](http://f4exb.free.fr/wsgc/proto_cpp/doc/html/index.html)
  * Essential functions _done (Dec 2012)_
    * Spent time fixing nasty bug with code using CUDA
    * Had to drop cublasIcamax method of finding maximum because norm works better than magnitude for weak signals
    * Piloted message correlator in CUDA _done_
  * Review unpiloted processing (for OOK modulation) _done (Jan 2013)_
    * Non-coherent detection (power) of OOK signal is much less performant. In the same conditions limit for all symbols decoded is S/N of -24dB for coherent BPSK and -6dB for non-coherent OOK
  * **Fading model finalization** _done (Jan 2013)_
    * See [Fading Models wiki page](http://code.google.com/p/wsgc/wiki/FadingModels) for some details
    * Watterson fading model works for any main sample rate fs in the limit of a Doppler spread of fs/250
  * Test with the fading model _started (Jan 2013)_
    * Conducted first informal tests. It appears that looking for a maximum on the cyclic averaging sum to determine symbol boundaries do not work very well. It does work when the signal is steady but fading introduces variation during the transmission of a symbol. Of course it depends on the fading rate but as the symbol period could be rather long (4s) it would be better to look for a specific synchronization scheme.
    * Hence more work is needed to design a proper synchronization scheme capable of establishing the epoch of the message and the symbols i.e. the point in time or in the samples continuous flow when the first PRN of the first symbol of the message is received.
    * **Synchronization scheme found** _(Jan 2013)_ It does produce errors at low S/N  and strong fading but in the same conditions the message would be unreadable anyway. So this seems to be the way to go. See [Synchronization Scheme wiki page](http://code.google.com/p/wsgc/wiki/SynchronizationScheme) for details.
      * Implement it in CUDA _basically works need debugging of samples processing in some cases_
      * Implement the external (to correlation) synchronization option in the CUDA message correlator _Done (Jan 2013)_
        * Bug identified see [Issue #1](https://code.google.com/p/wsgc/issues/detail?id=1). Implies a major rewrite in piloted message correlator CUDA flavour. _Issue closed (Jan 2013)_
    * Low-pass filtering: before jumping head on into a series of tests and benchmarks it could be useful to shape the signal closer to what it will be in the real world. Indeed we will want the signal to be low pass filtered using a raised cosine FIR filter for example. _Done implemented as an option (Jan 2013)_: a quick test shows it doesn't degrade performance.
    * Review and improve decoding box for synchronized messages _Done (Feb 2013)_
    * Review and improve strategy for the mitigation of correlation peak jitter due to multipath.
  * DBPSK study to mitigate the phase instability effect on long PRN sequences without resorting to completely incoherent OOK _Complete decoding of OOK and DBPSK implemented (Feb 2013)_
    * Quick results: we do observe frequency shift insensitivity. S/N limit is much worse than with BPSK but better than with OOK. In the same conditions (f<sub>s</sub>=4096Hz,  f<sub>c</sub>=1023Hz, m=10) limit is -13 dB instead of -24 dB.
  * Source coding a la JT65
    * Original scheme _done (Jul 2013)_
    * Extension to an optimal source coding for APRS
  * HD43 (by DJ5HG see Dubus 4/2009) is an interesting idea. What about a WSGC-43?
  * Forward Error Correction review and implementation (Reed-Solomon, convolutional, both?) _Started March 2013_
    * **I am spending quite a bit of brain energy on this topic and the project looks stalled although I am still quite active on it** Understanding how Reed-Solomon error correction works is quite involved in complex mathematics especially when aiming at advanced algorithms for soft decision decoding. In fact a proper complete decoding process should include everything and it is difficult to get relevant results and comparisons to the MFSK based technique without building a complete engine especially when soft-decision is involved. This implies understanding how the advandced soft-decision algorithms work. Exploiting just the results without going into the details of mathematical demonstrations is already quite complex. My intuition is that soft-decision should be involved when dealing with correlation peak or power spectrum density maxima. I expect it to make a significant difference as already illustrated with JT65 performance comparison between hard-decision and soft-decision Koetter-Vardy RS alogrithm.
    * _May 2013_ A separate project for the Reed-Solomon soft decoding (and Reed-Solomon encoding but this is more trivial) has been created as the [rssoft project](https://code.google.com/p/rssoft/) in Google Code. This is to make it usable outside the WSGC context. Therefore it will be in the form of a library. I think I have now understood how this all works at least functionally and I am confident to get something out!.
    * _June 2013_ Work on rssoft library is going on. A Galois Field layer has been written.
    * _July 2013_ **rssoft library is ready.** JT65 source encoding and soft decision Reed-Solomon options are implemented in wsgc\_test and wsgc\_generator prototypes ([see -j and -O options](https://code.google.com/p/wsgc/wiki/CppPrototype#Source_coding_options)). For now this is only tested for MFSK scheme (a la WSJT). Correlation and demodulation classes have to be adapted to soft-decision. However things are made simpler because the point is just to fill the reliability matrix with the correlation or demodulation values. The core of the soft-decision process is handled at the end of the reception period in a specific "decision box".
    * _August 2013_ A few benchmarks are carried out to exercise RSSoft in the MFSK domain that makes it easy to compare with WSJT/JT65. See [MFSK, JT65 and RSSoft study](https://code.google.com/p/wsgc/#MFSK,_JT65_and_RSSoft_study) paragraph above.
  * **_August 2013_ Due to the benchmarks that were run the correlation of Gold Code sequences just do not stand the comparison with plain MFSK when dealing with weak signal transmission. Therefore WSGC as such is dropped as a weak-signal transmission mode. No need to do more complete tests (see the top of this page for details)**

## Live software system ##

May be used for something else than WSGC.

  * Radio interface server
  * Computation engine server
  * Monitor and control client

## Radio hardware system ##

![http://ed06.ref-union.org/ACTIVITES/Dsc_7810.jpg](http://ed06.ref-union.org/ACTIVITES/Dsc_7810.jpg)

_No, it's not me on the picture. The radio gear will not be as massive but still..._

May be used for something else than WSGC.

  * Local oscillator
    * The local oscillator is very critical to the system and should achieve outstanding accuracy and stability
  * Receiver front end
    * Possibly around Analog Devices [ADRF6806](http://www.analog.com/en/rfif-components/modulatorsdemodulators/adrf6806/products/product.html) for 50 to 432 MHz bands
  * Transmitter front end
    * Possibly around Analog Devices [ADRF6755](http://www.analog.com/en/rfif-components/modulatorsdemodulators/adrf6755/products/product.html) for 144 MHz to 2.4 GHz bands

# About me #

http://ed06.ref-union.org/ACTIVITES/Doublier_200710.JPG

_F4EXB/P in action. This is an old picture, now I sit in the front of a computer..._

I am a licensed Amateur Radio operator with callsign F4EXB since July 2005. As many fellow "hams" I had interest in the radio business since a very young age. My grandfather who I didn't have the chance to meet was among the pioneers and participated to the first radio transmission between Panthéon and Tour Eiffel in Paris in the last years of the 19th century. I live in the sunny Riviera in the south of France. My QRA locator is JN33NN. My main interests in radio are weak signals and higher bands, above 144 MHz and up to 10 GHz with prospects to run 24, 47 and 76 GHz rigs. I am also interested in Laser communications especially in the infrared domain around 1064 and 1550nm where both clear terrestrial propagation windows and fairly accessible technology exists. With WSGC I might take interest in the lower HF bands (I have 80m band in mind) where it would be easier to make the first live tests.

My day job is in the computer science business.