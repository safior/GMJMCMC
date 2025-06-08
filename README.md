# GMJMCMC with correlation features

The folder "R" contains all the code used to run the base algorithm, and the extensions. Every file ending in ts, re and re2 are files created by me, and not a part of the base GMJMCMC package. 

As of now, there is quite a lot of duplicated code. Some functions within the files with ts, re or re2 endings are completely new, and some are modified. The extent of the modifications varies between simply including more input parameters, to substantial altercations.  

If this code was to be included in the GMJMCMC package, it should be cleaned up, and duplications should be removed. But to clearly indicate what I have done, and what I have not done, no files from the base GMJMCMC package have been removed, or edited in any substantial way. This should make it clear which functions I have created, and the extent to which they have been modified compared to the functions contained in the base GMJMCMC package; one can simply compare the function with ts, re or re2 ending with its counterpart without any of the three endings.
