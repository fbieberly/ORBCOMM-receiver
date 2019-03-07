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

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf - No longer available online. PDF is in the literature folder.  
http://www.decodesystems.com/orbcomm.html  
 

## Data format

In the data folder is a couple files of samples that I have recorded.  

The files are .mat files. They can be opened with MATLAB or Python (using SciPy's loadmat function).  

The files include metadata:  
- fc: center frequency  
- fs: sample rate  
- sats: a list of the names of the satellites overhead  
- tles: a list of lists of the tle lines for each satellite (in the order of the sats list)  
- timestamp: unix time of the start of the recording  
- samples: a numpy complex64 array of the samples  

