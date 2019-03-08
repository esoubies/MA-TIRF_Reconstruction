# MA-TIRF_Reconstruction

## Description

Multi-angle total internal reflection fluorescence microscopy (MA-TIRF) reconstruction algorithm which jointly performs volume reconstruction, deconvolution, and background estimation. More details can be found in the following papers:

<a href="https://www.nature.com/articles/s41598-018-36119-3" target="_blank">Nanometric Axial Resolution of Fibronectin Assembly Units Achieved with an Efficient Reconstruction Approach for Multi-Angle-TIRF Microscopy.</a>
Scientific Reports, 9-1, 2019.  <br />
Emmanuel Soubies, Agata Radwanska, Dominique Grall, Laure Blanc-Féraud, Ellen Van Obberghen-Schilling, and Sébastien Schaub.

<a href="https://hal.inria.fr/hal-02017862" target="_blank">Improving 3D MA-TIRF Reconstruction with Deconvolution and Background Estimation. </a>
Proc. ISBI, 2019. <br />
Emmanuel Soubies, Laure Blanc-Féraud, Sébastien Schaub, Ellen Van Obberghen-Schilling.

<p align="center">
<img src="https://github.com/esoubies/MA-TIRF_Reconstruction/blob/master/Images/data.png"/>
<img src="https://github.com/esoubies/MA-TIRF_Reconstruction/blob/master/Images/recons.png"/>
</p>

## Requirements

The code requires the GlobalBioIm library v1.1.1 (or more recent releases)  https://biomedical-imaging-group.github.io/GlobalBioIm/

## Repository content
* main function **RecTIRF.m** 
* function **TIRF_matrix.m** which creates the MA-TIRF matrix given physical parameters
* function **OptiADMMtirf.m** : ADMM class derived from the OptiADMM class of GlobalBioIm
* function **getColorCodedDepthFig.m** which generates the color-coded depth representation of the reconstructed image
* folder **Example** containing a simple example of use in **Script.m** 
