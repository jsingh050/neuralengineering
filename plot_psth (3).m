function [psth,tnew]=plot_psth(t,raster,tbin)
% Usage ... [psth,tnew]=plot_psth(t,raster,tbin)
%
% Generates a PSTH plot from Raster data based
% on time bin width tbin

if nargin==2,
  tbin=raster;
  raster=t;
  t=[1:size(raster,1)];
elseif nargin==1,
  error('  you must select a time-bin');
end;

dt=t(2)-t(1);

tnew=[t(1):tbin:t(end)];
for mm=1:length(tnew),
  search_i=find((t>=tnew(mm))&(t<(tnew(mm)+tbin)));
  psth(mm,:)=sum(raster(search_i,:)>0.5,1);
  %psth(mm,:)=sum(raster(search_i,:)>0.5,1)/tbin;
end;

plot(tnew,psth)
axis('tight'), grid('on'),
ylabel('PSTH (#/bin)'),

