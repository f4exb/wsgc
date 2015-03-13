# JT Modes galore! #

Exploring even more modes based on JT protocol



## MFSK variants for JT65 ##

![https://wsgc.googlecode.com/git/img/MFSK_JT65B_vs_MFSK_variants.png](https://wsgc.googlecode.com/git/img/MFSK_JT65B_vs_MFSK_variants.png)

This figure shows the decoding performance with variants stemming from original JT65B. MFSK description follows the wsgc\_test -d MFSK option format with the log2 of symbol bandwidth as first subparameter and the log2 of the symbol length as the second subparameter:
  * <font color='#000080'><b>JT65B</b></font>: This is the JT65B reference line. In fact it is MFSK:0,1 shifted by -0.5 dB as it fits nearly perfectly the documented JT65B line
  * <font color='#ff0000'><b>MFSK:1,0</b></font>: in fact this is the clone of JT65B in the wsgc\_test simulator but using RSSoft soft-decoding algorithm and a 4096/11025 frequency scaling factor from JT65B or a 11025/4096 time scaling factor as explained in the [KV vs RSSoft wiki page](https://code.google.com/p/wsgc/wiki/RSSoft_vs_KV).
  * <font color='#008000'><b>MFSK:1,1</b></font>: The symbol length has been increased from 1s to 2s
  * <font color='#C0C000'><b>MFSK:0,0</b></font>: The symbol bandwidth has been reduced from 2 Hz to 1Hz
  * <font color='#008080'><b>JT129</b></font>: Variant of the source coder that packs the 72 bits + 10 bits of the digits 5 and 6 of the locator (see next for the 10 bits) into 12 7-bit symbols with 2 bits spare. It uses an alphabet of 128 symbols and a RS(127,12) code instead of the RS(63,12). The codeword sent is 127 symbols long.
  * <font color='#800080'><b>JT257</b></font>: Variant of the source coder that packs the 72 bits into 9 8-bit symbols therefore using an alphabet of 256 symbols and a RS(255,9) code instead of the RS(63,12). The codeword sent is 255 symbols long.

We can draw the following conclusions:
  * We do not benefit of reducing symbol bandwidth as much as with increasing symbol length. MFSK:0,0 has a 1 dB gain over MFSK:1,0 but MFSK:1,1 has a 2 dB gain over MFSK:1,0
  * Using a longer RS(255,9) code instead of the original RS(63,12) is clearly not worth the extra message length. It is much more beneficial to increase symbol time

## Packing more bits ##

It would be nice to pack a few more bits in a JT65 message in particular to be able to transmit 6 digit locators. The two extra digits take 24\*24 = 576 positions so they need 10 bits to be encoded. Two extra 6-bit symbols could be added to the original JT65 protocol with 2 spare bits.

We saw in the preceding paragraph that increasing the n/k ratio of the RS(n,k) code does not yield tremendous improvement. Therefore we can expect that reducing this ratio will not hinder the performance too much. Provided we use a RS(63,14) code instead of the original RS(63,12) code we will be able to convey these two extra bits with no increase in codeword length.

Let's see how the RS(63.14) of the 6 digit locator mode that we call JT65\_L6 with MFSK:1,0 modulation compares with the original JT65 with MFSK:1,0:

![https://wsgc.googlecode.com/git/img/MFSK_JT65B_vs_JT65_L6.png](https://wsgc.googlecode.com/git/img/MFSK_JT65B_vs_JT65_L6.png)

JT65\_L6 (yellow curve) falls more abruptly than MFSK:1,0 (red curve) but is still very close at the falloff point of -22 dB.

Of course both curves are below the JT65B reference (blue curve) with the new JT65\_L6 being even further away. It could still be interesting to try this mode if 6 digit locators are really required.

## Here comes JT33 ##

We saw in the first paragraph that changes on RS length and ratio have a moderate influence on the overall performance compared to the increase in symbol length. The trouble with say doubling the symbol length is that it doubles the length of JT messages as well and therefore the time to complete a QSO.

Here the compromise is to shorten the code by using 5-bit symbols and a RS(31,k) code accordingly while doubling the symbol time. The result is that we send 31 symbols twice the length so it makes a 62 symbol period compared to 63. The resulting sent message is even a little shorter than JT65.

As there are 32 symbols we will call it JT33 in a similar way as JT65.

We still have to determine the value of k to fit our JT65 72 bits source encoded message. The equation is fairly easy: 15\*5 = 75. So we'll have 15 5-bit symbols with 3 spare bits. Therefore a RS(31,15) code is used.

While we are at it let's consider 6 digit locators. This adds 10 bits so exactly 2 symbols and we would then have 17 symbols and a RS(31,17) code would be used.

Below is the figure with the results:

![https://wsgc.googlecode.com/git/img/MFSK_JT65_vs_JT33.png](https://wsgc.googlecode.com/git/img/MFSK_JT65_vs_JT33.png)

Well this is what I call a good compromise that doesn't do any compromise on performance!

We see that JT33 (yellow curve) is even a tiny bit better than original JT65B (blue curve). JT33\_L6 (green curve) does a tiny bit worse but is clearly ahead of MFSK:1,0 (red curve) which is the original JT65B but with RSSoft decoder.

At least this is what it looks like without Doppler spread. Unfortunately it is very sensitive to Doppler spread and performance degrades in much larger proportions than with JT65 RS(63,12) classical scheme when Doppler fading is introduced.