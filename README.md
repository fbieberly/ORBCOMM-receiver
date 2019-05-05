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
[NumPy]: https://www.numpy.org/
[SciPy]: https://www.scipy.org/
[pyrtlsdr]: https://github.com/roger-/pyrtlsdr




## Getting started
  
### Offline recording and decoding  
1. First run the update_orbcomm_tle.py script to get the latest two-line elements for the orbcomm satellites.  
2. Update latitude and longitude of your receiver in record_orbcomm.py  
3. Record IQ data by running record_orbcomm.py  
4. Run file_decoder.py to decode a single recording file (defaults to the first file in the /data folder)  
  
  
### Real-time recording and decoding  
Not implemented yet.  
  
  
  
## References

I used these two resources as my primary references.  

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf
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

Look at the mat_file_explorer.py script to see an example of how to access the metadata.  

