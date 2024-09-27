function plot_raster(t,y)
% Usage ... plot_raster(t,y)
%
% Generates a stacked line plot of the columns of y
% vs time t. If no time vector is supplied, indeces are used

if nargin==1,
  y=t;
  t=[1:size(y,1)];
end;

for mm=1:size(y,2),
  if mm==1, hold('on'), end;
  plot([t(1) t(end)],mm+[0 0],'k-'),
  plot(t,mm+single(y(:,mm)>0.5)*0.5,'k-'),
end;
hold('off'),

ylabel('Raster'),

if nargin==1,
  xlabel('Index #'),
else,
  xlabel('Time'),
end;

axis('tight'), grid('on'),

