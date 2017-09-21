function [ X, exitFlag ] = alignbyRegion( X, axisscale, ppm_max, ppm_min, type)
%     Align spectra using cross-correlation in a specific region
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     axisscale: spectral axisscale in ppm
%
%     OUTPUT
%     X: ppm referenced spectra
%     exitFlag: exception code

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

if isempty(ppm_max) || isempty(ppm_min),
    ppm_max=axisscale(1);
    ppm_min=axisscale(end);
end
    
if ppm_max<=ppm_min,
    exitFlag=1;
    disp('Maximum ppm must be higher than the minimum ppm');
    return;
else
    [~,index_min]=min(abs(axisscale-ppm_max));
    [~,index_max]=min(abs(axisscale-ppm_min));
end

[r,c] = size(X);
x=zeros(r, c);

wb = waitbar(0,'Aligning. Please wait...');
waitbar(0/r);

for i=1:r,
    if i==1,
        %Use the median to align the first sample
        ref = median(abs(real(X(:, index_min:index_max))),1);
    else
        ref = abs(median(real(x(1:i, index_min:index_max))));
    end
    %substract offsets
    ref=ref-min(ref);
    ref=ref';
    sample = abs(real(X(i, index_min:index_max))); 
    sample=sample-min(sample);    
    sample=sample';
    aux = [ref, sample];    
    [c, l]=xcorr(aux);
    [~, ix]=max(c(:,2));
    if strcmp('all', type),
        x(i,:)=circshift(X(i, :),[0 l(ix)]); 
    else
        x(i,:)=[X(i, 1:index_min-1)  circshift(X(i, index_min:index_max),[0 l(ix)]) X(i, index_max+1:end)];
    end
    waitbar(i/r-1);
end
X=x;
close(wb) 
end

