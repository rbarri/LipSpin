function [ x ] = changeSpectralRange( x, ranges, type, dim )
%     Resize spectra according to indicated regions
%
%     INPUT
%     x: struct with (at least) the following fields:
%        data: [NxM] spectra, where N is the number of spectra and M the spectral data points
%        axisscale: [1xN] cell with spectral axisscales
%     ranges: [1x2] matrix of N included spectral ranges in the form [max_ppm, min_ppm] 
%     type:
%          'i'- make resized spectra including only selected ranges
%          'e'- make resized spectra excluding selected ranges
%     dim: dimension

%     OUTPUT
%     x: resized x

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

if (ranges(1)<ranges(2))
    disp('Set ranges from higher to lower ppms')
    return;
end
l=size(x.axisscale{dim},2);
interval=[];
for i=0:length(ranges)/2-1;
    [~,low_index]=min(abs(x.axisscale{dim}-ranges(i*2+1)));
    [~,high_index]=min(abs(x.axisscale{dim}-ranges(i*2+2)));
    interval=[interval low_index:high_index];        
end
if type == 'e',    
    interval= delsamps([1:l]', interval)';
end
ndim=ndims(x.data);
if dim==2,
    if ndim==3,
        x.data=x.data(:,interval,:);
    else
        x.data=x.data(:,interval);
    end
elseif dim==3,
    x.data=x.data(:,:,interval);
end
x.axisscale{dim}=x.axisscale{dim}(interval);
end

