function [ FWHH, skew] = getFWHH( X, axisscale )
%     Return the full width at half height and the swekness of a single peak in X
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     axisscale: spectral axisscale in ppm
%
%     OUTPUT
%     FWHH: full width at half height
%     skew: skewness as the difference in simmetry from center at HH

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

%get the real part
X=real(X);
nbSamples=size(X, 1);
FWHH=zeros(nbSamples, 1);
skew=zeros(nbSamples, 1);

if license('test','signal_toolbox'),
    for i=1:nbSamples,
        [intensities, indexes]=findpeaks(X(i,:), 'SORTSTR', 'descend'); 
        center=axisscale(indexes(1));
        hw=intensities(1)/2;
        [~, hw1]=min(abs(X(i,1:indexes(1))-hw));
        hw1=axisscale(hw1);
        [~, hw2]=min(abs(X(i,indexes(1)+1:end)-hw));
        hw2=axisscale(hw2+indexes(1)+1);        
        FWHH(i)=(abs(hw1-center)+abs(hw2-center));
        skew(i)=max(abs(hw1-center), abs(hw2-center))/min(abs(hw1-center), abs(hw2-center));
    end
else
    for i=1:nbSamples,
        aux=diff(X(i,:),2);
        [~, index]=min(aux);
        index=index+1;
        center=axisscale(index);
        intensity=X(i,index);
        hw=intensity/2;
        [~, hw1]=min(abs(X(i,1:index)-hw));
        hw1=axisscale(hw1);
        [~, hw2]=min(abs(X(i,index+1:end)-hw));
        hw2=axisscale(hw2+index+1);        
        FWHH(i)=(abs(hw1-center)+abs(hw2-center));
        skew(i)=max(abs(hw1-center), abs(hw2-center))/min(abs(hw1-center), abs(hw2-center));
    end
end
end

