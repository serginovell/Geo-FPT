# GEO-FPT
Geometrical Fitted Bispectrum model based on Perturbation theory. Presented in the paper . 

The C code for the GEO-FPT bispectrum model is found in C_lib/B02_rsd.c. 

First, one should compile the Makefile in ./C_lib, just typing "make" in a terminal located in that folder. If there are any errors, check for the installed packages or the directories in the Makefile.

Once that is done, one can run it from Python, using the wrapper function bk_multip from the interface.py file.

The needed input quantities for bk_multip are:

tr, tr2,tr3,tr4: Triangle configurations to be computed for each multipole: (0), (200), (020), (002) respectively. They have to be in an array of size (Ntriangles,3)

kp: Input power spectrum k bins
pk: Input power spectrum. GEO-FPT uses the non-linear matter power spectrum as input, which is computed first with PTcool and then combined with BRASS (this code provides the neede functions from BRASS to combine the power spectrum, in the Pdd_P02.c file. All credit for this is for Héctor Gil-Marín.

cosm_par: An array with the cosmological parameters, in the following order: [ $\sigma8$ , $f$, $\alpha_\parallel$, $\alpha_\bot$, $b_1$ , $b_2$ , $A_P$, $\sigma_P$, $A_B$, $\sigma_B$]. See the paper in [] for the details on each parameter.

The redshift at which we want to compute the bispectrum

fit_full: Set it to 1 in order to use the model fitted with the full data-vector B02. Otherwise, it will use the model calibrated only with B0.

If instead one wants to use the function directly in C, one could use the void ext_bk_mp function in B02_rsd.c. The input parameters are very similar, with the exception that in C the power spectrum and its k-bins have to be inputted as their log10. Also, the number of elements of the power spectrum have to be specified in kp_dim, as well as the number of triangles of the multipoles (0), (200), (020), (002), in the quantities num_tr, num_tr2,num_tr3,num_tr4 respectively.

The most relevant data-vectors used in the paper [] can be found in the folder ./input_data/. The Quijote and Nseries bispectrum triangles, means, covariances and means+shot-noise are found in the file quijote_dvs.npz: labelled as "quij_tr...", "quij_mean...","quij_cov...","quij_sn...", with z indicating the redshift. As for the triangles, tr0 being the monopole and 200 quadrupole, and tr020 being the triangles for the quadrupole 020. The _sn data-vectors are the bispectrum without subtracting the shot-noise. So in order to obtain the shot-noise, one could do quij_sn-quij_mean.

All is analogous for Nseries.
