Simple example:
```julia
filename = "/Users/ferreol/Data/SPHERE/PSF/right_cropped_reduced_SPHER.2021-11-24T03:32:33.190_IRDIS_OBJECT,FLUX_0.837464s_5f_DP_0_BB_Ks_POLARIMETRY,CORONOGRAPHY.fits";

hdu = read(FitsHeader, filename)
img = read(FitsFile, filename; ext=1)
w = read(FitsFile, filename; ext=2)
hdr = read(FitsHeader, filename)



Diameter =8.2 # diametre des UT en m
RAD2MAS = 1 / π *180 *60*60*1000 
λ = 2e-6 #(doit être fixée d'apres les keywords)
pixscal=hdr["PIXSCAL"].value()

pupil_radius = pixscal / RAD2MAS  * Diameter / λ
badpix = w .==0;
BoucheTrou!(img,badpix,pupil_radius)
```