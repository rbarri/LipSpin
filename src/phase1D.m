function [ X, ph0, ph1 ] = phase1D( X, axisscale, mode, ranges, resolution, iterations, auto)
%   Batch correction of phase distortions in 1D-NMR spectra
%
%   INPUT
%   X: [NxM] 1D spectra, where N is the number of spectra and M the spectral data points
%   axisscale: spectral axisscale in ppm
%   mode: 
%      0 (only 180 phase is checked based on maximum absorbtive part in the selected ranges)
%      1 (minimize entropy using Chen Li method)(Journal of Magnetic Resonance 158 (2002) 164–168)
%      2 (minimize ls error respect to a flat baseline)
%      3 (force ranges to be symmetric) (not implemented yet)
%   ranges: [Nx2] matrix of N included spectral ranges in the form [max_ppm, min_ppm] 
%   resolution: length (in ppm) for data binning
%   iterations: maximum number of iterations before stopping
%   auto: automatic range selection (0,1) (only for mode 2)
%
%   OUTPUT
%   X: baseline-corrected spectra
%   ph0: angle correction for ph0
%   ph1: angle correction for ph1

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

nbSamples=size(X, 1);

if isreal(X),
    X=hilbert(X);
end

[Xaux newAxisscale]=bindata(X, 0.02,  axisscale);

ph0=zeros(1,nbSamples);
ph1=zeros(1,nbSamples);

interval=[];

if ~isempty(ranges),
    for i=1:size(ranges,1),
        %remove wrong ranges
        if ranges(i,1)<=ranges(i,2),
            ranges(1,:)=[];
        else
            interval=[interval find(newAxisscale<ranges(i,1) & newAxisscale>(ranges(i,2)))];
        end
    end
end

subranges=[];
if isempty(ranges)
    interval=1:length(newAxisscale);
    if auto==1 && mode==2,
       nBins=round((newAxisscale(1)-newAxisscale(end))/resolution);
       nBins=round((interval(end)-interval(1))/nBins);
       subranges=interval(1):nBins:interval(end);
    end    
end

wb = waitbar(0,'Phasing. Please wait...');
waitbar(0/nbSamples);

switch mode,
    case 0 %Check if 180º correction is neccesary
        invert=-sum(real(Xaux(:, interval)), 2)>sum(real(Xaux(:, interval)), 2);
        for i=1:nbSamples,
            if invert(i),
                ph0(i)=pi;
                ph1(i)=0;
                X(i,:) = setPhase(X(i,:), pi, 0);
            end
        end
    case 1 %Apply phase correction based on entropy
        for i=1:nbSamples
            best=realmax;
            bestph0=realmax;
            bestph1=realmax;
            j=1;
            iNoImprovement=0;
            while j<=iterations,                
                disp(['Phasing by entropy. Iteration number: ' num2str(j)]);               
                ph0=(rand(1)*2*pi);
                ph1=(rand(1)*2*pi);   
                [fac f]=fminsearch(@entropyPhase, [ph0 ph1], optimset('TolX',1e-4), Xaux(i,:), interval);
                %delimit phases between 0 and 2pi to avoid 2pi multiples in
                %order to increment iNoImprovement
                ph0Aux=(fac(1)-floor(fac(1)/(2*pi))*2*pi); 
                ph1Aux=(fac(2)-floor(fac(2)/(2*pi))*2*pi);                
                diff_ph0=abs(bestph0-ph0Aux);
                diff_ph1=abs(bestph1-ph1Aux);
                if best>f && (diff_ph0>1e-3 || diff_ph1>1e-3),
                    best=f;
                    bestph0=ph0Aux;
                    bestph1=ph1Aux;
                    iNoImprovement=0;
                else
                    iNoImprovement=iNoImprovement+1;
                end
                if iNoImprovement==3,
                    disp(['Phasing by entropy. No more improvements after ' num2str(iNoImprovement) ' iterations']);                
                    break;
                end
                j=j+1;
            end
            waitbar(i/nbSamples);
            %set phases in the [-pi pi] range
            if bestph0>pi,
                bestph0=bestph0-2*pi;
            end
            if bestph1>pi,
                bestph1=bestph1-2*pi;
            end
            X(i,:) = setPhase(X(i,:), bestph0, bestph1);
        end
    case 2 %Apply phase correction based on flattening the baseline      
        lb=[-pi -pi -1]; %Lower bound (-180º)
        ub=[pi pi 1]; %Upper bound (180º)
        x0=[0 0 0]; %Initial solution   
        for i=1:nbSamples,
            %First, apply automatic range selection
            if auto==1,
                if ~isempty(subranges)
                   for j=1:length(subranges)-1,
                       Xdiff=diff(real(Xaux(i, subranges(j):subranges(j+1))));
                       stdX(j)=std(Xdiff);
                   end
                   p60stdX=prctile(stdX, 60);
                   intervalAux=interval;
                   for j=1:length(subranges)-1,
                       %less than 40% of the spectral regions is expected to
                       %have peaks (could not work in all cases)
                       if stdX(j) > p60stdX;
                           interval(interval>=subranges(j) & interval<=subranges(j+1))=[];
                       end
                   end
                   pulseOnes=ones(1,nBins*2);
                   specBreaks=double(diff(interval)>1);
                   specBreaks=conv(specBreaks, pulseOnes);
                   if any(specBreaks),
                        interval=interval(logical(specBreaks(nBins:end-nBins+1)));
                   else
                       %Do nothing
                   end
                   if isempty(interval)
                       interval=intervalAux;
                   end
                end
            end
            
            ref=zeros(1, length(interval));         
            xaux=[];
            xaux(1,:)=Xaux(i,interval)/max(abs(Xaux(i,interval)));
            xaux(2,:)=interval;
            xaux(3,1)=length(Xaux(i,:));
            best=realmax;
            bestph0=realmax;
            bestph1=realmax;
            j=1;
            iNoImprovement=0;            
            while j<=iterations,
                disp(['Phasing by flattening the baseline. Iteration number: ' num2str(j)]); 
                x0(1)=(rand(1)*2*pi)-pi;
                x0(2)=(rand(1)*2*pi)-pi;               
                x0(3)=0;
                options=optimset('TolX',1e-6, 'display', 'off');
                [fac res] = lsqcurvefit(@flatBaselinePhase, x0, xaux, ref, lb, ub, options); 
                diff_ph0=abs(bestph0-fac(1));
                diff_ph1=abs(bestph1-fac(2));
                if best>res && (diff_ph0>1e-3 || diff_ph1>1e-3),
                    best=res;
                    bestph0=fac(1);
                    bestph1=fac(2);
                    iNoImprovement=0;
                else
                    iNoImprovement=iNoImprovement+1;
                end
                if iNoImprovement==3,
                    disp(['Phasing by flattening the baseline. No more improvements after ' num2str(iNoImprovement) ' iterations']);                
                    break;
                end
                j=j+1;       
            end
            waitbar(i/nbSamples);
            X(i,:) = setPhase(X(i,:), bestph0, bestph1);
        end
    case 3 
end

close(wb) 
end

function f = entropyPhase(x, X, interval)
    ph0=x(1);
    ph1=x(2);     
    %dephase 
    X = setPhase(X, ph0, ph1);
    X=real(X(interval));
    order=1;
    x_der=abs(diff(X,order));
    p1=x_der/sum(x_der);
    %Calculation of Entropy
%     lp1=length(p1);
    p1(1,(p1(1,:)==0))=1;
%     for i=1:lp1
%         if (p1(1,i)==0)%in case of ln(0)
%             p1(1,i)=1; 
%         end
%     end
    h1  = -p1.*log(p1);
    f  = sum(h1);        
end

function f = flatBaselinePhase(x, X)
    ph0=x(1);
    ph1=x(2);
    offset=x(3);
    indexes=X(2,:);
    N=X(3,1);
    f=setPhase(X(1,:), ph0, ph1, N, indexes);
    f=real(f)+offset;
    f((isnan(f)))=0;
end