function [Xout, newAxisscale] = bindata(X, binSize, axisscale)
%     Spectral binning  
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     binSize: bin size in ppm 
%     axisscale: non-binned spectral axisscale in ppm
%
%     OUTPUT
%     Xout: binned spectra
%     newAxisscale: binned spectral axisscale in ppm

%     Copyright (C) 2017 Rubén Barrilero Regadera
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

%First, set the limits to integrate the bin
xbinLimits=axisscale(1):-binSize:axisscale(end);
newAxisscale=zeros(1,length(xbinLimits)-1);

for i=1:length(xbinLimits)-1,
    [~,ppmmax]=min(abs(xbinLimits(i)-axisscale));
    [~,ppmmin]=min(abs(xbinLimits(i+1)-axisscale));
    newAxisscale(i) = mean([axisscale(ppmmax) axisscale(ppmmin)]);
end

if ndims(X) == 3,
    Xout=zeros(size(X,1), length(newAxisscale), size(X,3));
    for j=1:size(X,3),
        for k=2:length(xbinLimits)
            Xout(:, k-1, j)=nanmean(squeeze(X(:, (axisscale>xbinLimits(k) & axisscale<=xbinLimits(k-1) ), j)), 2);
        end
    end    
else
    Xout=zeros(size(X,1), length(newAxisscale));
    for k=2:length(xbinLimits)
        Xout(:, k-1)=nanmean(X(:, (axisscale>xbinLimits(k) & axisscale<=xbinLimits(k-1) )),2);
    end
end
end
