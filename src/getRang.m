function [ interval ] = getRang( axisscale, maxppm, minppm )
%     This function returns the indexes of the interval between maxppm and
%     minppm included in axisscale
%
%     INPUT
%     axisscale: spectral axisscale in ppm
%     maxppm: upper chemical shift for aligned region
%     minppm: lower chemical shift for aligned region
%
%     OUTPUT
%     interval: indexes of axisscale between maxppm and minppm

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

[~,max_index]=min(abs(axisscale-minppm));
[~,min_index]=min(abs(axisscale-maxppm));

interval=min_index:max_index;
end

