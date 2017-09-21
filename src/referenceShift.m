function [X, shift] = referenceShift(X, axisscale, ppm_ref, shift_tolerance)
%     Move spectra to match maximum peak signal in a spectral range with ppm_ref
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     axisscale: spectral axisscale in ppm
%     ppm_ref: chemical shift where the maximum signal of a reference peak has to be moved 
%     shift_tolerance: left and rigth boundaries from ppm_ref within the peak is sought
%
%     OUTPUT
%     X: aligned spectra
%     shift: shift (in units) 

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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>


if ndims(X)<3,
    X_=X;
else
    X_=sum(X,3);
end
    
    
[r,c]=size(X_);

shift=zeros(r,1);

[~, minindex]=min(abs(axisscale-(ppm_ref+shift_tolerance)));
[~, maxindex]=min(abs(axisscale-(ppm_ref-shift_tolerance)));

wb = waitbar(0,'Referencing. Please wait...');
waitbar(0/r);

for i=1:r,
    bjump=0;
    for j=1:(maxindex-minindex+1), 
        %find the number of equal points 
        repet=numel(find(real(X_(i, minindex:maxindex))==real(X_(i, minindex+j))));
        %if one of the repeated points is the maximum, select it directly
        %without findpeaks
        if repet>=2 && real(X_(i, minindex+j))==max(real(X_(i, minindex:maxindex)))
            bjump=1;
            locs(1)=j;
            break;
        end
    end
    %Find reference index
    [~, index]=min(abs(axisscale-ppm_ref));    
    if bjump==0,
        %find peaks with default
        [~,locs] = findpeaks(abs(real(X_(i, minindex:maxindex))), 'SORTSTR', 'descend');
    end
    %Find first peak.
    if isempty(locs),
        %Sometimes, the maximum of a peak is shared by two or more peaks or
        %no peak is found, in that case choose maximum (not necesarely the expected solution)
        disp(['peaks not found for sample: ' num2str(i)]);
        [~, locs(1)]=max(abs(real(X_(i, minindex:maxindex))));
    end
    shift(i,1)= index-(locs(1)+minindex-1);   
    X(i, :, :)=circshift(X(i, :, :), [0 shift(i,1) 0]);
    waitbar(i/r);
end

close(wb)
