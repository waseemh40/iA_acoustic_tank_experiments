Description of algorithims:

avergae_envelope -> for tag and combining 8 pulses. Uses echo rejection and excludes frequencies. Fixed experiment time of 3 minutes.
avergae_envelope_v2 -> same as above except number of sets is removed and scans file till end. Takes more time for calculations.

comparison_filter -> used for tag, single pulse approach. Fixed experiment timem of 3 minutes. Echoi rejection and frequency exclusion schemes are same as avergae_enevelope.
comparison_filter_v2 -> same as above except that it takes time as argument.
comparison_filter_v3 -> same as v2 but eho rejection time is changed to 800ms. This is trimmed for transducer NOT tag. Also number of sets check removed in here. Also in next commit (in 3rd I guess), a plot for freq shifd/speed will also be added in this version.


function_generator -> Slices time axis in Fs time and does fft. Nothing threshold related things. 

Rest are earlier version used in development.




