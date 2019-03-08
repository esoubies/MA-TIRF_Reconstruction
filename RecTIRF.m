%--------------------------------------------------------------------
% This script performs 3D reconstrcution from MA-TIRF measurements
%    - Data-term: Least-Squares 
%    - regul: TV or Hessian-Schatten norm
%
% References:
% [1] Emmanuel Soubies, Agata Radwanska, Dominique Grall, Laure Blanc-Feraud, Ellen Van Obberghen-Schilling, 
%     and Sebastien Schaub. "Nanometric axial resolution of fibronectin assembly units using multi-angle 
%     total internal reflection fluorescence microscopy", Scientific Reports, 9-1, 2019. 
% [2] Emmanuel Soubies, Laure Blanc-Feraud, Sebastien Schaub and Ellen Van Obberghen-Schilling, 
%     "Improving 3D MA-TIRF Reconstruction with Deconvolution and
%     Background Estimation", Proc ISBI 2019
%
% This script makes use of the GlobalBioIm library (v 1.1.1 or more recent releases) :
%     https://biomedical-imaging-group.github.io/GlobalBioIm/
%
% Before to run it set the following papareters:
% -- General 
%  pathToScript=   -> Path to main script TIRFrecons.m (the present one)
%  dataname=       -> File name to stack TIRF acquisitions (.mat file)
%  psfName=        -> File name to PSF (.mat file). Only used if deconv is activated
%
% -- MA-TIRF reconstruction
% - One can directly provide the TIRF matrix 
%  tirfOpname=     -> File name to TIRF matrix (nb angles x nb z, .mat file)
% - Used only if tirfOpname=[] (empty)
%  zmax=           -> Max z value for reconstruction (micron)
%  rz=             -> Desired axial resolution  (micron)
%  rxy=            -> Lateral resolution  (micron)
%  ni=             -> refractive index of the incident medium 
%  nt=             -> refractive index of the transmitted medium
%  lamb=           -> exitation wavelength (micron)
%  Na=             -> numerical aperture
% - Other parameters
%  deconv=         -> Boolean true to include joint deconvolution, see [1,2]
%  valback=        -> Background value that will be subtracted to data
%  estiback=       -> Boolean true to include background estimation, see [2]
%  lambBack=       -> Regularization parameter for the background estimation, see [2]
%
% -- ADMM Parameters 
%  mu=             -> Regularization parameter (can be an array to loop)
%  muL1=           -> Regularization parameter for L1 prior, see [2](if zero, no L1 prior is used)
%  maxIt=          -> Max iterations
%  Reg=            -> Choice regul: 1 for TV, 2 for Hessian-Schatten 
%  rhoDT=          -> rho parameter (ADMM) associated to data term 
%  rhoReg=         -> rho parameter (ADMM) associated to the regularization
%  rhoPos=         -> rho parameter (ADMM) associated to non-negativity
%--------------------------------------------------------------------

%% Copiright (2018) Emmanuel Soubies (esoubies@gmail.com)
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------
global isGPU

%% Reading data
tmp=load(dataname);
fld=fieldnames(tmp);
y=getfield(tmp,fld{1});
y=y-valback;maxy=max(y(:));y=y/maxy;
szData=size(y);

%% Build PSF and Illumination matrix
% -- PSF
if deconv
    tmp=load(psfName);
    fld=fieldnames(tmp);
    psf=getfield(tmp,fld{1});
end

% -- Illumination matrix
if ~isempty(tirfOpname)
    tmp=load(tirfOpname);
    fld=fieldnames(tmp);
    A=getfield(tmp,fld{1});
    A=A/max(A(:));
    norColA=ones(1,size(A,2));    % Matrix given by the user : no column normalization
    assert(size(y,3)==size(A,1),'The number of slices in acquisition must be equal to the number of rows of the TIRF matrix');
else
    tmp=load(anglName);
    fld=fieldnames(tmp);
    angl=getfield(tmp,fld{1});
    z=[rz/2:rz:zmax]';dz=ones(size(z))*rz;
    A=TIRF_matrix(z,dz,angl,lamb,ni,nt);
    norColA=sqrt(sum(A.^2,1));
    A=A./repmat(norColA,[size(A,1) 1]);
end
szRecons=[szData(1:2),size(A,2)];

%% Conversion CPU/GPU is necessary
if deconv, psf=gpuCpuConverter(psf); end
y=gpuCpuConverter(y);
A=gpuCpuConverter(A);

%% Common operators and costs
% -- Data term
TIRF=LinOpBroadcastMatrix(A,szRecons,3);
L2=CostL2([],y)*TIRF;
if deconv
    H=LinOpConv(repmat(fftn(psf),[1 1 szRecons(3)]),1,[1 2]);        % Convolution operator
else
    H=LinOpIdentity(szRecons);
end
% -- Regularization
if Reg==1
    Opreg=LinOpGrad(szRecons,[1,2]);               % TV regularizer: Gradient operator
    Freg=CostMixNorm21(Opreg.sizeout,4);           % TV regularizer: Mixed Norm 2-1
elseif Reg==2
    Opreg=LinOpHess(szRecons,'circular',[1 2]);    % Hessian-Shatten: Hessian Operator
    Freg=CostMixNormSchatt1(Opreg.sizeout,1);      % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
end
% -- Non-Negativity constraint
Id=LinOpIdentity(szRecons);
if muL1~=0
    pos=CostNonNeg(szRecons)+muL1*CostL1(szRecons);
else
    pos=CostNonNeg(szRecons);
end

%% TIRF Reconstruction 
for ii=1:length(mu)
    % ADMM
    FF=[{L2},{mu*Freg},{pos}];
    HH=[{H},{Opreg},{Id}];
    rho=[rhoDT,rhoReg,rhoPos];
    Opt=OptiADMMtirf([],FF,HH,rho,estiback,lambBack);
    if muL1~=0, costPos=[1 2 4]; else costPos=[1 2]; end;
    Opt.OutOp=OutputOpti(1,[],round(maxIt/10),costPos);
    Opt.ItUpOut=round(maxIt/10);      % call OutputOpti update every ItUpOut iterations
    Opt.maxiter=maxIt;                % max number of iterations
    Opt.run(zeros_(szRecons));         % run the algorithm
    
    % -- Display color-coded image
    if isGPU==1
        xopt=gather(Opt.xopt)./repmat(reshape(norColA,[1,1,length(norColA)]),[size(Opt.xopt,1),size(Opt.xopt,2),1]);
        if estiback
            back=gather(Opt.back)/sqrt(szData(3));
        end
    else   
        xopt=Opt.xopt./repmat(reshape(norColA,[1,1,length(norColA)]),[size(Opt.xopt,1),size(Opt.xopt,2),1]);  
        if estiback
            back=Opt.back/sqrt(szData(3));
        end
    end
    getColorCodedDepthFig(xopt,rz);
    if estiback
        figure; imagesc(back); axis image; axis off; 
        title('Estimated Background');set(gca,'FontSize',16);
        colormap gray; colorbar;
    end
end