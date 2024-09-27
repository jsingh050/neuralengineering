function y=imshift2(x,y0,x0,zero_flag)
% Usage ... y=imshift2(im,x0,y0)
%
% 2D shift in fourier domain

verbose_flag=0;

if nargin<4, zero_flag=0; end;
  
% shift the function by the transit time in Frq Domain
xfft = fftshift(fft2(x));
xphase = [-size(x,2)/2:size(x,2)/2-1]*2*pi/size(x,2);
xphase = ones(size(x,1),1)*xphase;
%show(xphase), pause,
xphase = x0*xphase;
xphase = exp(-i.*xphase );
yphase = [-size(x,1)/2:size(x,1)/2-1]*2*pi/size(x,1);
yphase = yphase'*ones(1,size(x,2));
%show(yphase), pause,
yphase = y0*yphase;
yphase = exp(-i.*yphase );
xfft2 = xfft.*xphase;
xfft2 = xfft2.*yphase;

% carry out the convolution by multiplying in frq. domain
y = real(ifft2(ifftshift(xfft2)));

if verbose_flag, disp(sprintf('  imshift = [%.4f, %.4f]',x0,y0)); end;

if zero_flag,
  if y0>1, 
    y(1:floor(y0),:)=0;
  elseif y0<-1,
    y(end+floor(y0)+1:end,:)=0;
  end;
  if x0>1, 
    y(:,1:floor(x0))=0;
  elseif x0<-1,
    y(:,end+floor(x0)+1:end)=0;
  end;
end;

 
if (nargout==0),
  show(y)
  clear y
end;


