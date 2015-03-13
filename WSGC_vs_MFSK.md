# WSGC vs MFSK #

This wiki page discusses the performance of WSGC correlation scheme with respect to simple MFSK. In all cases except the JT65B reference a Reed-Solomon code is used with the RSSoft soft-decoder. JT65B uses KVASD and is 0.5 dB better as discussed [here](https://code.google.com/p/wsgc/wiki/RSSoft_vs_KV).

The first results are disappointing to say the least. WSGC is way behind MFSK so it's clearly not worth the effort.

For now only AWGN has been used to simulate noise, I am waiting for more comparisons in adverse Doppler fading conditions before dropping this project. I have not much hope that it will do any better than MFSK.

**This is it, when it comes to weak signals nothing beats just plain non-coherent MFSK**

Below the sad news:

The MFSK figures are taken from [JT modes galore!](https://code.google.com/p/wsgc/wiki/MFSK_JT_modes_galore)

![https://wsgc.googlecode.com/git/img/WSGC_vs_MFSK.png](https://wsgc.googlecode.com/git/img/WSGC_vs_MFSK.png)


WSGC is the light blue line. It is using 1023 chip PRNs, 1s per PRN hence 4s per symbol. In JT65 terms the symbol time would be 4 times longer at 1,488s instead of 0.372s.

The figure speaks by itself. No need for further comments...