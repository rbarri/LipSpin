function [ noiseLevel dynRange] = getNoiseLevel (X)
%     This function returns a level of noise that is calculated from the
%     standard deviation of the noisy spectral region in X after baseline
%     supression.
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral
%     data points
%
%     OUTPUT
%     noiseLevel: RMS noise level
%     dynRange: difference between maximum and minimum points in X

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

%get noise region length
N=size(X,2);
baseline=[];

if license('test','signal_toolbox'),
    windowSize=round(N*0.3);
    if ~rem(windowSize,2),
        windowSize=windowSize+1;
    end
    baseline=sgolayfilt(X,2,windowSize);
elseif license('test', 'Curve_Fitting_Toolbox'),
    baseline=smooth(X, round(N*0.3), 'rloess');
    baseline=baseline';    
else
    %use a randomly sampled subset to calculate the interpolation
    idx = randperm(N);
    idx(round(N*0.3):end)=[];
    idx=sort(idx);
    xfirst=[];
    xlast=[];
    x=X(1,idx);
    %Not necessary but add estimation of border points
    if idx(1)~=1,
        xfirst=mean(X(1, 1:idx(1)));
        idx=[1 idx];
    end
    if idx(end)~=N,
        xlast=mean(X(1, idx(end):end));
        idx=[idx N];
    end
    x=[xfirst x xlast];
    x=(x-mean(X))/std(X);
    p=polyfit(idx, x, 3);
    baseline = polyval(p,1:N);
    baseline=(baseline.*std(X))+mean(X);
end
dynRange=max(max(X))-min(min(X));
X=X-baseline;
noiseLevel=std(X);
end