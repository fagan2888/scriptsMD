function [danom,dmean,coeff] = removeharm(mvarced,sampfreq)
%function [danom,dmean,coeff] = removeharm(mvarced,sampfreq)
%  remove first 3 annual harmonics (stored in dmean) from input signal.
%   [danom,dmean] = removeharm(mvarced)
%  input signal (MVARCED) must be a column vector!
%  SAMPFREQ is sampling interval (e.g. 1/5 for once every 5 days)

sz=size(mvarced);
ntimes=sz(1);

year=[1:ntimes]./(365*sampfreq);
D=mvarced;
plotctr='n';

% Plot a longitude-time Hovmuller diagram
if plotctr=='y'
    lon = linspace(0, 360, 512);
    lat = linspace(-90, 90, 256);
    figure;
%     set(gcf,'pos',[655   172   358   510])
    clf
    lati = 256/2 + 10;
%     if lat(lati)<0
%         H = ['^oS'];
%     else
%         H = ['^oN'];
%     end
    contourf(lon,year,squeeze(mvarced(:,lati,:)));
%     shading flat;
    colorbar;
%     set(gca,'tickdir','out')
%     titlestr = ['longitude and time at latitude ' ...
%         num2str(abs(lat(lati))) H];
%     title(titlestr)
end


%  Specify design matrix:
%  Remove first 3 harmonics:
Eann = [ones(size(year(:))) cos(2*pi*year(:)) sin(2*pi*year(:)) ...
    cos(2*pi*2*year(:)) sin(2*pi*2*year(:)) ...
    cos(2*pi*3*year(:)) sin(2*pi*3*year(:))];
Dann = zeros(size(D));
for m=1:sz(2)
    d = D(:,m);
    coeff = Eann\d(:);
    dfit = Eann*coeff;
    Dann(:,m) = dfit;
end
%   dmean = reshape(Dann,[nlats nlons ntimes]);
dmean=Dann;
%   D = D-Dann;
%   danom = reshape(D,[nlats nlons ntimes]);
danom=D - Dann;


if plotctr=='y'
    % replot Hovmuller
    figure
    contourf(lon,year,squeeze(danom(:,lati,:)));
%     shading flat;
    colorbar;
%     set(gca,'tickdir','out')
%     title(titlestr)
%     disp(' ');disp('Re-plotted Hovmuller with 3 harmonics removed');
%     figure(gcf)
end
