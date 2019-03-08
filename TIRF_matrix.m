function I=TIRF_matrix(z,dz,alpha,lamb,ni,nt)
%--------------------------------------------
% function I=TIRF_matrix(z,dz,alpha,lamb,ni,nt)
%
% Computes the TIRF matrix for a given setting
%
% Inputs : z     -> z-discretization (vector containing the center of each bin)
%          dz    -> associated size of the bins (vector same size as z)
%          alpha -> vector containing the used incident angles (radian)
%          lamb  -> exitation wavelength
%          ni    -> refractive index of the incident medium (in general glass 1.518)
%          nt    -> refractive index of the transmitted medium (in general water 1.333)
%
% Output : I     -> TIRF matrix defined by :
%                        I(l,j) = I0(l)/p(l) [ exp[-(z(j)+dz(j)/2)*p(l)] - exp[-(z(j)-dz(j)/2)*p(l)] ]
%                   where the angles are indexed by l and the z by j. I0 denotes the intensity at the 
%                   interface and p is the inverse of the penetration depth.
%
% Copyright (2015) Emmanuel Soubies (esoubies@gmail.com)
%--------------------------------------------

% -- Get the number of angles and the lenght of z
nb_a=length(alpha);
nb_z=length(z);

% -- Check the constitency of the length of vector z and dz
if (length(z)~=length(dz))
	error('In TIRF_matrix : inputs z and dz must have the same length !');
end

% -- Declaration of the result matrix
I=zeros(nb_a,nb_z);

% -- Initialization critical angle and the inverse of the penetration depth
ac=asin(nt/ni); % critical angle
p=(4*pi*ni)./lamb*sqrt(sin(alpha).^2-(nt/ni)^2);

% -- I0 Axelrod
% Polarisation p
pol_p=4*cos(alpha).^2.*(2*sin(alpha).^2-(nt/ni)^2)./((nt/ni)^4*cos(alpha).^2+sin(alpha).^2-(nt/ni)^2);
% Polarisation s
pol_s=4*cos(alpha).^2./(1-(nt/ni)^2);
I0=0.25*pol_p+0.75*pol_s;

zl=z-dz/2; % left bound of the bins (in the z-discretization)
zr=z+dz/2; % right bound of the bins (in the z-discretization)
if (~isempty(find(zl<0)) || ~isempty(find(zr<0)))
	error('In TIRF_matrix : incoherent values of z and dz !');
end

for l=1:nb_a
	I(l,:)=I0(l)/p(l)*[exp(-p(l)*zl)-exp(-p(l)*zr)];
end

end
