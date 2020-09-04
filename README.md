# Death adder lure reflectance

This repository contains the code used to plot diet data and reflectance measurements in southern death adders ("Acanthophis antarcticus").

R-scripts were run on a local machine where the reflectance measurements were stored as text files output from OceanView Software. Measurements were taken on an Ocean Optics spectrometer (MAYA2000 Pro).
reflectance measurement was expressed relative to a Spectralon 99% white reflectance standard (WS-1-SL, Ocean Optics, Dunedin, Fl, USA) and the probe (QR200-7-UV-BX 200µm) was mounted to maintain a constant angle (90°) and distance (4 mm) from the surface of the skin. 

The 'data' directory contains a sample info (.csv) with information on taxanomy, snout-vent legnth, mass, and tail lure morphology. 
This directory also contains a directories with reflectance measurements (.txt), and diet information (.csv) collated from the literature.

The 'plots' directory contains the resulting figures output by the Rscripts.