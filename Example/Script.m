%-----------------------------------------------------------
% File containing parameters for MA-TIRF reconstruction
%     see script TIRFscript.m
%
% Copyright (2018) Emmanuel Soubies (esoubies@gmail.com)
%-----------------------------------------------------------
clear;%close all

% -- To run on GPU (0: CPU / 1: Matlab Parrallel Computing Toolbox / 2: CudaMat) 
%    Note: when using Matlab Parrallel Computing Toolbox (useGPU(1)), resolution with HessianSchatten Norm 
%          may not be optimal because svd computations are done with a mex function (see Cost/CostUtils/HessianSchatten) 
%          on the cpu (so lost of time in transfert CPU/GPU).
useGPU(0)

%% General 
pathToScript='../';                   % Path to main script TIRFrecons.m
dataname='Example/Acquisitions';      % File name stack TIRF acquisitions (.mat file)
psfName='Example/psf';                % File name to PSF (.mat file). Only used if deconv is activated

%% MA-TIRF reconstruction
% - One can directly provide the TIRF matrix 
tirfOpname=[];               % File name to TIRF matrix (nb angles x nb z, .mat file)
% - Used only if tirfOpname=[] (empty)
anglName='Example/angles';   % File to angles (radian, .mat file)
zmax=0.300;                  % Max z value for reconstruction (micron)
rz=0.02;                     % Desired axial resolution  (micron)
rxy=0.106;                   % Lateral resolution  (micron)
ni=1.518;                    % refractive index of the incident medium 
nt=1.34;                     % refractive index of the transmitted medium
lamb=0.491;                  % exitation wavelength (micron)
Na=1.33;                     % numerical aperture 
% - Other parameters
deconv=1;                    % Boolean true to include joint deconvolution
valback=200;                 % Background value that will be subtracted to data
estiback=0;                  % Boolean true to include background estimation
lambBack=1e2;                % Regularization parameter for the background estimation

%% ADMM
mu=[5e-4];         % Regularization parameter (can be an array to loop)
muL1=0;            % Regularization parameter for L1 prior (if zero, no L1 prior is used)
maxIt=50;          % Max iterations
Reg=1;             % Choice regul: 1 for TV, 2 for Hessian-Schatten 
rhoDT=1e-1;        % rho parameter (ADMM) associated to data term 
rhoReg=1e-1;       % rho parameter (ADMM) associated to the regularized
rhoPos=1e-1;       % rho parameter (ADMM) associated to non-negativity

%% Run reconstruction algo
run([pathToScript,'RecTIRF']);
