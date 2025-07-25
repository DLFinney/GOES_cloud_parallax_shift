# GOES_cloud_parallax_shift
Code which adjusts GOES cloud fields for parallax shift, and regrids them.

Developed as part of the DCMEX project https://cloudsense.ac.uk/dcmex/

These are being shared as inspiration for others coding similar problems, and for others to check what has been 
used for my work. The scripts weren't developed specifically for sharing so you may find they are not simply plug and play
but I've done my best to make the steps clear.

__Some notes on parallax shift__

In a satellite image, the apparent position of a cloud can have an error of several kilometres if the parallax effect is not accounted for. This effect arises because, when calculating the latitude and longitude coordinates of pixels in a satellite image, an assumption must be made for the height of that point above the centre of the earth. A simple assumption of mean sea level can introduce errors on the order of kilometres. Errors are largest for higher clouds towards the edge of the satellite image. Since DCMEX looks at deep convective clouds, and because New Mexico is towards the edge of the GOES 16 satellite image, the parallax error could not be missed when trying to relate aircraft tracks to cloud images from the satellite. Since DCMEX studies growing congestus with highly variable cloud heights, an adaptive approach for parallax correction was developed. This method individually corrected each pixel in each timestep by utlising the GOES cloud-top height product.


__Other python code that may be suitable for the task__

SatPy https://satpy.readthedocs.io/en/stable/api/satpy.modifiers.parallax.html

__Other python code that may be suitable for the task__
This methodology has been developed and used as part of the DCMEX project.
Finney et al. (2025) Microphysical fingerprints in anvil cloud albedo. ACP. https://doi.org/10.5194/egusphere-2025-1227 
Finney et al. (2024) Deep Convective Microphysics Experiment (DCMEX) coordinated aircraft and ground observations: microphysics, aerosol, and dynamics during cumulonimbus development. ESSD. https://doi.org/10.5194/essd-16-2141-2024
