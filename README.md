# wavelet-opt
Matlab Files to optimize a wavelet function using GA and lifting schemes (LS), which are an implementation of the WT using filter banks. The LS have two FIR filter, U and P, that correspond to the wavelet and scalling functions in the traditional WT. The GA optimize U and P based on a given cost function, for example, class separability, mutual information, etc. This method has been presented among other references in [1,2].

[1]R D Pinzon-Morales and Yutaka Hirata. 2012. Customization of
wavelet function for pupil fluctuation analysis to evaluate levels of 
sleepiness. In Proceedings of the 11th international conference on Teleco
mmunications and Informatics, Proceedings of the 11th international 
conference on Signal Processing (SITE'12),115-120. 

[2]R D Pinzon Morales et al. Novel signal-dependent filter bank method for
identification of multiple basal ganglia nuclei in Parkinsonian patients
2011 J. Neural Eng. 8 036026
