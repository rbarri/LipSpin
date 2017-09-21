function [ f ] = getApodizationFunction(type, timeScale, lb, gm, ssb)
%     Create common apodization functions
%
%     INPUT
%     type: apodization function
%           'lor' - lorentzian function 
%           'gau' - gaussian function 
%           'sin' - shifted sine-bell function
%     timeScale: axisscale for FID
%     lb: exponential broadening factor used for lorentzian and gaussian windows
%     gm: gaussian broadening factor factor in gaussian windows
%     ssb: order in sine function
%
%     OUTPUT
%     f: apodization function

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

nPoints=length(timeScale);
switch type
    case 'lor'
        f=exp(-(timeScale)*lb*pi);
    case 'gau'
        a=pi*lb;
        b=-(a/(2*gm*timeScale(end)));
        f=exp(-a*(timeScale)-(b*(timeScale.^2)));
    case 'sin'
        if ssb<=1,
            f=sin((pi()-0)*((0:nPoints-1)/(nPoints-1)));
        else
            ssb=floor(ssb);
            f=sin((pi()-pi()/ssb)*((0:nPoints-1)/(nPoints-1))+pi()/ssb);
        end
end
end

