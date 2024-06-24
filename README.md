# AMSR2-CSAT-Combined
Optimal estimation combined algorithm for AMSR2-CloudSat coincident overpass data. 

This is the source code for the optimal estimation (OE) combined algorithm for drizzle retrievals used in the study Jones and Kummerow 2024 (manuscript submitted for peer review). Source code for retrieval is in Fortran90 and jupyter notebooks contain code for creating the figures in the paper. Figures can also be downloaded in their native jpeg format for reproduction purposes with proper citations.

Some important notes:
1. Due to data limitations, many of the lookup tables are omitted, but the source code is possible to run without the tables using some of the native subroutines included.
2. Other subroutines or data not included were not part of the work and were either created for other works or are publicly available via other sources. Those are:
       1) Water vapor EOFs were stored in binary format and were created by David Duncan for Duncan and Kummerow (2016).
       2) Gaseous absorption lookup tables were created using MonoRTM. Subroutine for reading in the tables is given in the code for the user to understand the format of the tables. MonoRTM is available from Atmospheric and Environmental Research (https://github.com/AER-RC).
       3) Gamma particle size distribution lookup table and snow scattering lookup table are not provided, but were created using the mie.f subroutine and the read routines for them are provided in dsd.f90 to aid the user in reconstructing them.
       4) Surface emissivity model used was FastEM6 (Liu et al. 2011).
       5) HDF5 libraries are not included.
3. SOFTWARE IS PROVIDED FOR REFERENCE PURPOSES ONLY AND COMES WITH NO WARRANTY, IMPLIED OR OTHERWISE, AND NO GUARANTEE OF CROSS-PLATFORM COMPATIBILITY.


Any questions should be directed to Spencer Jones: spencer.jones@colostate.edu


References:


Duncan, D.I. and C.D. Kummerow, 2016: A 1DVAR retrieval applied to GMI: Algorithm description, validation, and sensitivities. Journal of Geophysical Research: Atmospheres. 121, 7415-7429. https://doi.org/10.1002/2016/JD024808

Liu, Q., Weng, F., and S.J. English, 2011: An Improved Fast Microwave Water Emissivity Model. IEEE Transactions on Geoscience and Remote Sensing, 49, 4, 1238-1250. https://doi.org/10.1109/TGRS.2010.2064779
