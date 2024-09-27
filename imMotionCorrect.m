function [xmc,ims]=imMotionCorrect(im1,im0,parms)
% Usage ... [xmc,ims]=imMotionCorrect(im1,im0,parms)
%
% Cross-correlation based approach for 2D rigid-body translation-only
% motion correction of images in a data set. If im1 and im0 are single
% images, im1 is corrected relative to im0. If im0 is a single number
% then the reference image is selected from the 3D set of images in im1.
% parms = [upsample norm_flag]
% when upsample is the upsampling scale factor to detect differences in
% motion, norm_flag is to cut values above and below the range of im1

im1=squeeze(im1);

im1sz=size(im1);
if length(im1sz)==2, im1sz(3)=1; end;

%im1sz,

if nargin<3, parms=[0.01 0]; end;
if length(parms)==1, parms(2)=0; end;

upf=parms(1);
zero_flag=parms(2);

if length(im0)==1, im0=im1(:,:,im0); end;

im0=im_smooth(im0,1);

im0_xc=fftshift(ifft2(fft2(im0).*conj(fft2(im0))));
[im0_xc_max_i, im0_xc_max_j]=find(im0_xc==max(im0_xc(:)));

[xx,yy]=meshgrid([1:size(im0_xc,1)], [1:size(im0_xc,2)]');

ims=zeros(im1sz);
for mm=1:im1sz(3);
  im1_orig=im1(:,:,mm);
  im1(:,:,mm)=im_smooth(im1(:,:,mm),1);
  im01_xc=fftshift(ifft2(fft2(im0).*conj(fft2(im1(:,:,mm)))));

  [im01_xc_max_i, im01_xc_max_j]=find(im01_xc==max(im01_xc(:)));
  i1_shift_est1 = [im0_xc_max_i-im01_xc_max_i  im0_xc_max_j-im01_xc_max_j];

  if 0,
  figure(gcf),
  subplot(211), plot([im0_xc(im01_xc_max_i,:)' im01_xc(im01_xc_max_i,:)']),
  axis tight, grid on,
  subplot(212), plot([im0_xc(:,im01_xc_max_j)  im01_xc(:,im01_xc_max_j) ]), 
  axis tight, grid on,
  drawnow, pause,
  end
  
  % let's try to use simple interp to find the sub-pixel amount
  [xxi,yyi]=meshgrid([-2:upf:2]+im01_xc_max_i, [-2:upf:2]'+im01_xc_max_j);
  im01_xci=interp2(xx,yy,im01_xc.',xxi,yyi,'cubic').';

  [im01_xci_max_i, im01_xci_max_j]=find(im01_xci==max(im01_xci(:)));
  if length(im01_xci_max_i)>1, im01_xci_max_i=im01_xci_max_i(1); end;
  if length(im01_xci_max_j)>1, im01_xci_max_j=im01_xci_max_j(1); end;
  
  i1_shift_est2 = [im0_xc_max_i-xxi(1,im01_xci_max_i)  im0_xc_max_j-yyi(im01_xci_max_j,1)];

  ims(:,:,mm) = imshift2(im1_orig, -i1_shift_est2(1), -i1_shift_est2(2));
  
  xmc(mm,:) = i1_shift_est2;
end;

if (zero_flag&(nargout==2)),
  ims(find(ims<ims_min))=ims_min;
  ims(find(ims>ims_max))=ims_max;
end
