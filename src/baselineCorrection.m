function [X] = baselineCorrection (X, axisscale, mode, ranges, resolution, orderPoly)
%     Remove baseline in spectra by interpolating functions to a set of spectral points in user-defined ranges  
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     axisscale: spectral axisscale in ppm
%     mode: subtraction/interpolation mode:
%           0-Median subtraction
%           1-Cubic spline data interpolation
%           2-Piecewise Cubic Hermite Interpolating Polynomial
%           3-N-order polynomial interpolation
%     ranges: [Nx2] matrix of N included spectral ranges in the form [max_ppm, min_ppm] 
%     resolution: spectral binning factor in ppm
%     orderPoly: order for polynomial interpolation
%
%     OUTPUT
%     X: baseline-corrected spectra

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

xPoints=[];
yPoints=[];

for i=1:size(ranges,1),
    %remove wrong ranges
    if ranges(i,1)<=ranges(i,2),
        ranges(1,:)=[];
    end
end

auto=0;
if isempty(ranges),
    %if ranges is empty take the length of the data
    auto=1;
    ranges(1,1)=axisscale(1);
    ranges(1,2)=axisscale(end);
end

wb = waitbar(0,'Correcting baseline. Please wait...');
waitbar(0/nbSamples);

xPoints=zeros(nbSamples,0);
for i=1:size(ranges,1),
    nBins=round(max((ranges(i,1)-ranges(i,2))/resolution, 2));    
    subranges=linspace(ranges(i,1), ranges(i,2), nBins);
    index=zeros(1,2);
    newPointIndex=size(xPoints, 2)+1;
    for j=1:nbSamples,
        ii=newPointIndex;
        [~,index(1)]=min(abs(axisscale-subranges(1)));  
        for k=1:length(subranges)-1,
            [~,index(2)]=min(abs(axisscale-subranges(k+1)));            
            stdBins(j, ii)=std(real(X(j, index(1):index(2))));
            yPoints(j, ii)=median(real(X(j, index(1):index(2))), 2);
            [~, xPoints(j, ii)]=min(abs(yPoints(j, ii)-real(X(j, index(1):index(2)))));
            xPoints(j, ii)=xPoints(j, ii)+index(1)-1;
            index(1)=index(2)+1;
            ii=ii+1;
        end
    end
end
nbPoints=size(xPoints,2);
[xPoints sortedIndexes]=sort(xPoints, 2);
iOffset=0:nbPoints:nbSamples*nbPoints;
iOffset=repmat(iOffset(1:end-1)', [1 nbPoints]);
sortedIndexes=sortedIndexes+iOffset;
if nbSamples>1,
    yPoints=yPoints';
end
yPoints=yPoints(sortedIndexes);

for i=1:nbSamples,
    isOK=0;
    iter=0;
    xNoise.data=X(i,:);
    xNoise.axisscale{2}=axisscale;
    xNoise=changeSpectralRange( xNoise, [axisscale(1) (axisscale(1)-0.5)], 'i', 2 );
    noiseLevel = getNoiseLevel(xNoise.data);
    if auto==1,
        %Regions with std below a threshold are kept
        selPoints=find(stdBins(i, :)<median(stdBins(i, :))+0.15*std(stdBins(i, :)));
        xPts=xPoints(i, selPoints);
        yPts=yPoints(i, selPoints);
        %Points with median below a threshold are kept
        selPoints=find(yPts<median(yPts)+std(yPts));
        xPts=xPts(selPoints);
        yPts=yPts(selPoints);
    else
        xPts=xPoints;
        yPts=yPoints;
    end
    while ~isOK && iter<5,
        %Several iterations are run to reduce negative values
        isOK=1;
        switch mode
            case 0
                baseline=median(yPts);        
            case 1
                baseline=spline(axisscale(xPts), yPts, axisscale);
            case 2
                baseline=pchip(axisscale(xPts), yPts, axisscale);              
            case 3
                p=polyfit(axisscale(xPts), yPts, orderPoly);
                baseline=polyval(p, axisscale);               
        end
        tempX=X(i,:)-(baseline+sqrt(-1)*baseline);
%         if mode~=0 &&  mode~=3,
%             newxPts=[];
%             newyPts=[];
%             for j=1:length(xPts)-1,
%                 [minValueY, minValueX]=min(real(tempX(xPts(j):xPts(j+1))));
%                 tempY=real(X(i,minValueX+xPts(j)));
%                 if minValueY<(-10*noiseLevel),
%                     isOK=0;
%                     newxPts=[newxPts minValueX+xPts(j)];
%                     newyPts=[newyPts tempY];
%                 end
%             end
%             xPts=[xPts newxPts];
%             yPts=[yPts newyPts];
%             [xPts, newOrder]=unique(xPts);
%             yPts=yPts(newOrder);
%             iter=iter+1;
%         end
    end
    waitbar(i/nbSamples);
    X(i,:)=tempX;
end
close(wb)
end