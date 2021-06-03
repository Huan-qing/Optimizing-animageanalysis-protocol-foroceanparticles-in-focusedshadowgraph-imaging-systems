# Image-preprocessing
MATLAB code for FoSI

# pipeline
1. Organize image data (delete_not_ctd.m)

   Remove the images not match with CTD

2. Schlieren filter (sch_filt.m)
   
   Filter out those not wanted schlieren images
   
3. Illumination correction (correction_function_noclip.m)

4. Subtraction (diff2pchip.m)

5. Particle analysis (diff2pchip.m)

6. Match the image data with CTD (diff2pchip.m)

7. Estimate data to the depth level (diff2pchip.m)
