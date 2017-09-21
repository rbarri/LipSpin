function X = setPhase(X, ph0, ph1, N, n)
%   Apply phase correction to a complex 1D spectrum
%  
%   INPUT:
%   X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%   ph0: zero-order phase correction
%   ph1: first-order phase correction
%   N: number data points of spectra
%   n: row vector with indexes of spectral points used for phasing.
%
%   OUTPUT:
%    X: phase-corrected spectrum

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

    if nargin<=3,
        N=size(X,2);
        n=1:N;
    end
    ph=ph0+ph1*(n/N);
    x = real(X).*cos(ph)-imag(X).*sin(ph);
    y = imag(X).*cos(ph)+real(X).*sin(ph);  
    X=x+1i.*y;
end