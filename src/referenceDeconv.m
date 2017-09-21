function [ X ] = referenceDeconv( X, axisscale, center, limits, mode, SF, refIndex, include )
%   Correct lineshape distortions using a reference spectra
%
%   INPUT:
%   X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%   axisscale: spectral axisscale in ppm
%   center: signal center
%   limits: signal limits [high PPM, low PPM]
%   mode: 1-Signal in dataset, 2-TMS signal
%   SF: frequency of spectrometer 
%   refIndex: sample index in X used as reference
%   include: column vector of included samples in X (0,1)
%   
%   OUTPUT:
%   X: lineshape-corrected spectra

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

%Get indexes for limits
[~, minIndex]=min(abs(axisscale-limits(1)));
[~, maxIndex]=min(abs(axisscale-limits(2)));
%Limits for baseline adjust
marginPoints=round((maxIndex-minIndex)/8);
minIndexPolyMargin=marginPoints;
maxIndexPolyMargin=(maxIndex-minIndex)-marginPoints;

switch mode
    case 1
        %Signal in dataset
        ref=real(X(refIndex,:));
        ref=ref(1, minIndex:maxIndex);
        %Adjust baseline at the extremes of the selected region
        refMargins=real(ref(1, [1:minIndexPolyMargin maxIndexPolyMargin:length(ref)]));
        f=polyfit(axisscale([minIndex:minIndex+minIndexPolyMargin minIndex+maxIndexPolyMargin:maxIndex]), refMargins, 1);
        baseline=polyval(f, axisscale(minIndex:maxIndex));
        refAux=ref-baseline;
        %Complete the original spectral length with zeros
        ref=zeros(1, size(X,2));
        ref(1, minIndex:maxIndex)=refAux; 
        %Reverse fourier transform
        ref=ifft((ref));
        ref = ref(1:end/2);
        ref(2:end)=2*ref(2:end);
    case 2
        %TMS signal
        omega=(center*SF/1e6 + abs(axisscale(end))*SF/1e6);
        axisscaleHz=(axisscale+abs(axisscale(end)))*SF/1e6;
        int=1;       
        delta=(limits(1)-limits(2))*SF/1e6;
        if delta>120 %include satellites
            omega(2)=omega(1)+59.25-0.00155*SF/1e6;
            omega(3)=omega(1)-59.25-0.00155*SF/1e6;
            int(2)=0.00555;
            int(3)=0.00555;
        end  
        omega(4)=omega(1)+3.31-0.00013*SF/1e6;
        omega(5)=omega(1)-3.31-0.00013*SF/1e6;
        omega(6)=omega(1)+1.37-0.0002*SF/1e6;
        omega(7)=omega(1)-1.37-0.0002*SF/1e6;
        omega(8)=omega(1)+1.16-0.0002*SF/1e6;
        omega(9)=omega(1)-1.16-0.0002*SF/1e6;
        omega(10)=omega(1)+0.95-0.0002*SF/1e6;
        omega(11)=omega(1)-0.95-0.0002*SF/1e6;
        omega(12)=omega(1)+0.74-0.0002*SF/1e6;
        omega(13)=omega(1)-0.74-0.0002*SF/1e6;

        int(4)=0.0247;
        int(5)=0.0247;
        int(6)=0.0016;
        int(7)=0.0016;
        int(8)=0.0047;
        int(9)=0.0047;
        int(10)=0.0047;
        int(11)=0.0047;
        int(12)=0.0016;
        int(13)=0.0016;  
        
        peakPars.hwhh=0.4;
        peakPars.gauss=0;
        refAux=zeros(1,length(axisscaleHz));
        for n=1:length(omega)
            peakPars.int=int(n);
            peakPars.pos=omega(n);
            if peakPars.int>0
                refAux=refAux+voigtProfile(peakPars, axisscaleHz);
            end
        end
        ref=zeros(1, size(X,2));
        ref(1, minIndex:maxIndex)=refAux(1, minIndex:maxIndex);         
        ref=ifft((ref));
        ref=ref(1:end/2);
        ref(2:end)=2*ref(2:end);       
end

for i=1:length(include),
    target_long=real(X(include(i), :));
    target_long_aux=target_long;
    target_short=target_long_aux(1, minIndex:maxIndex); 
    targetMargins=real(target_short(1, [1:minIndexPolyMargin maxIndexPolyMargin:length(target_short)]));
    %Adjust baseline at the extremes of the selected region
    f=polyfit(axisscale([minIndex:minIndex+minIndexPolyMargin minIndex+maxIndexPolyMargin:maxIndex]), targetMargins, 1);
    baseline=polyval(f, axisscale(minIndex:maxIndex));
    target_short_aux=target_short-(baseline);    
    %Complete to the original spectral length with zeros
    target_short=zeros(1, size(X,2));
    target_short(1,minIndex:maxIndex)=target_short_aux;             
    target_short=ifft(target_short);
    target_short = target_short(1:end/2);
    target_short(2:end) = 2*target_short(2:end);     
    target_short=target_short./ref;
    %normalization
    target_short=target_short/target_short(1);    
    target_long_aux=ifft(target_long);
    target_long=target_long_aux;
    target_long = target_long(1:end/2);
    target_long(2:end) = 2*target_long(2:end);     
    %Apply the quotient to a limited region of the FID in
    %order to reduce wiggles
    invOffset=1;
%     weights=tukeywin(length(target_long)*2, 0.95)';
    weights=parzenwin(length(target_long)*2)';
    weights=weights(end/2+1:end);
    inv_weights=invOffset-weights;
    %Uncomment this line and comment the following to avoid weighting 
%     target_long=(target_long./target_short);
    target_long=(target_long./target_short).*weights + target_long.*inv_weights;   
    target_long(1)=target_long(1)*0.5;
    %Fourier transform      
    X(include(i), :)=fft(target_long, length(axisscale));
end
end

