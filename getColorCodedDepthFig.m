function [proj,colBar]=getColorCodedDepthFig(im,dz,prt)
%------------------------------------------------------
% function [proj,colBar]=getColorCodedDepthFig(im,dz,prt)
%
% Generates a color-coded depth representation of the given image
% (im) with the given axial resolution (dz). The parameter prt is a 
% bollean (default true) to display the result.
% This function returns both the color-coded depth representation (proj)
% as well as the colorbar (colBar).
%
% Copyright (2018) Emmanuel Soubies (esoubies@gmail.com)
%------------------------------------------------------

if nargin <3
    prt=1;
end

sz=size(im);
map=jet(sz(3));
map=map(end:-1:1,:);
R=zeros(sz);G=zeros(sz);B=zeros(sz);
for i=1:sz(3)
    R(:,:,i)=im(:,:,i)*map(i,1);
    G(:,:,i)=im(:,:,i)*map(i,2);
    B(:,:,i)=im(:,:,i)*map(i,3);
end
proj=cat(3,cat(3,sum(R,3),sum(G,3)),sum(B,3))/sz(3);
proj=proj/max(proj(:));

tmp=reshape(jet(sz(1)),sz(1),1,3);
colBar=repmat(tmp(end:-1:1,:,:),1,20,1);
if prt
    figure;
    image([colBar,ones([sz(1),10,3]),proj]);axis image; %axis off;
    ylabel(['0 - ',num2str(sz(3)*dz)]);
    set(gca,'Fontsize',16);set(gca,'xtick',[]);set(gca,'ytick',[])
end
end