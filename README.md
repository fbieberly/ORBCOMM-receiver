# orbcomm_decoder
A software receiver for ORBCOMM satellite transmissions.
  

  
## Description

This is a software receiver for decoding packets from ORBCOMM satellites. 
While all packets can be decoded and logged, only the Ephemeris packet is 
decoded and the results can be plotted.  

I am writing this decoder as an instructional personal project. Hopefully it
can be used by others to learn about designing satellite communication 
receivers.  
If you want a more full-featured ORBCOMM receiver please check out:  
https://www.coaa.co.uk/orbcommplotter.htm  
http://f6cte.free.fr/index_anglais.htm  



## Dependencies

Written for py2.7.  

I use [pyrtlsdr] to record the RF signal with an RTLSDR receiver.  
[NumPy] and [SciPy] are used for signal processing.  
[PyEphem] is used to calculate Az/El and doppler shift of the satellites.  

pip install pyrtlsdr, numpy, scipy, pyephem  



[PyEphem]: https://rhodesmill.org/pyephem/index.html
[NumPy]: www.numpy.org/
[SciPy]: https://www.scipy.org/
[pyrtlsdr]: https://github.com/roger-/pyrtlsdr


## References

I used these two resources as my primary references.  

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf  
http://www.decodesystems.com/orbcomm.html  
 
