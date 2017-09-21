function [vP] = voigtProfile(peakPars, axisscale)
%   Generate voigt profile (Mixed Lorentzian-Gaussian)
%
%   INPUT
%   peakPars: structure 4 element vector with peak parameters:
%      int: intensity
%      pos: center in ppm
%      hwhh: half width at half height in ppm
%      gauss: guassian fraction (0-1).
%   axisscale: spectral axisscale in ppm
%
%   OUTPUT
%   vP: resulting voigt profile

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

lor= 1./(pi()*peakPars.hwhh*(1+((axisscale - peakPars.pos)/peakPars.hwhh).^2)); %lorentzian shape 
gau=exp(-(((axisscale - peakPars.pos)/peakPars.hwhh).^2)/(2/(sqrt(2*log(2))^2))); %gaussian function, where sigma in simple equation is replaced by the HWHH of the voigt profile. It keeps the FWHH indicated no matter the G/L ratio
gau   = peakPars.gauss*gau;
lor   = (1-peakPars.gauss)*lor;
vP   = peakPars.int*(gau + lor); %weighted combination
%Scale to math the sought intensity
vP=vP/max(vP); 
vP=vP*peakPars.int;
end
