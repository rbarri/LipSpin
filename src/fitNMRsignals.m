function [ output_args ] = fitNMRsignals( X, axisscale, noiseLevel, mets, freq, BGtype, constraint_factor, maxError, tries, plotFinal)
%     least squares curve fitting of 1D NMR using multiplet models with voigt functions or spectral templates 
%
%     INPUT
%     X: [NxM] spectra, where N is the number of spectra and M the spectral data points
%     axisscale: spectral axisscale in ppm
%     noiseLevel: estimation of RMS noise level
%     mets: structure of substructures containing parameters for individual signals. Substructures have to be named as 'metx' and 'specx' for voigt and spectral models, respectively, with x ranging from 1 to the maximum number of signal elements.
%           Common elements for 'metx' and 'specx':
%                pos: signal position (in ppm)
%                fwhh: full width at half signal height (in Hz)
%                gauss: Gaussian fraction (0-1)
%                protons: number of H's providing that signal
%           Specific elements for 'metx':
%                jcop: J-coupling constant (in Hz)
%                mult: signal multiplicity (1,2,3...)
%           Specific elements for 'specx':
%                xaxis: axisscale of spectral template
%                data: spectral data points
%                stdPos: signal position in spectral template
%                hlim: upper chemical shift limit for spectral template (ppm)
%                llim: lower chemical shift limit for spectral template (ppm)
%     freq: spectrometer frequency in Hz
%     BGtype: bitwise variable, each bit activates one type of background: 
%             0x01 = 1st degree polynomial
%             0x02 = Gaussian
%             0x04 = Cosine background with first degree polynomial
%     constraint_factor: row vector containing constrain factors for met parameters  ranging from 0 to inf (0 means less flexibility) 
%                        Vector elements: 1-Intensity (not used), 2-position, 3-hwhh, 4-Jcoupling, 5-Gaussian
%     maxError: stop criterion for fitting based on maximum error allowed
%     tries: maximum number of fitting runs without reaching stop criteria
%     plotFinal: display fitting solution in a program window
%
%     OUTPUT
%     output_args: structure with the following fields:
%                  flag: exit condition for lsqcurvefit
%                  F: estimated spectrum  
%                  s: optimised signal parameters
%                  error: fitting RMSE
%                  h: figure handler
%                  pureSpectra: estimated spectra for individual signals
%                  baselineSpectra: estimated baseline

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


%select real part of the spectrum
X=real(X);

%Calculate the total number of variables in order to customize the fitting
submets=fieldnames(mets);
nbnVariables=length(submets)*5;

%Maximum intensity normalization to make intensity comparable with other
%fitting variables
scale_factor_intensities=max(abs(X));
X=X/scale_factor_intensities;
noiseLevel=noiseLevel/scale_factor_intensities;
scale_factor=0;

%Set the minimum change for partial derivatives in the non-linear least squares method
DiffMinChange=1e-4;
MinConstraintRange=DiffMinChange+DiffMinChange/100;

%% Generic boundaries for metabolites
%Set boundaries for intensity
[~,dynRange]=getNoiseLevel(X);
if dynRange>10*noiseLevel
    lbI=noiseLevel*5;
else
    lbI=noiseLevel;
end
ubI=1;
if ubI<lbI,
    ubI=1;
    lbI=0.001; %Noise is as higher as signal
end

%% Generic boundaries for background signals
%Set boundaries for intensity
lbIbckg=0;
ubIbckg=10;

%Set boundaries for chemical shift
lbppmbckg=-1;
ubppmbckg=1;

%Set boundaries for width
lfwhhbckg=(15*1e6)/freq;
ufwhhbckg=(200*1e6)/freq;

%Set boundaries for gaussian percentage
ugaussFactorbckg=1;
lgaussFactorbckg=0;

%% Initialize fitting parameters
[a,b]=size(constraint_factor);
if a==1 || b==1,
    constraint_factor=repmat(constraint_factor,1,length(submets));
end
for i=1:length(submets),
    fields=getfield(mets, submets{i});
    if strncmp('bckg', submets{i}, 4)
        s0(1,i)=lbIbckg + (ubIbckg-lbIbckg)*rand(1);
        ub(1,i)=ubIbckg;
        lb(1,i)=0;
        s0(2,i)=fields.pos;
        ub(2,i)=ubppmbckg+fields.pos;
        lb(2,i)=lbppmbckg+fields.pos;
        s0(3,i)=(fields.fwhh*1e6)/freq;
        ub(3,i)=ufwhhbckg;
        lb(3,i)=lfwhhbckg;
        ub(4,i)=ugaussFactorbckg;
        lb(4,i)=lgaussFactorbckg;        
        s0(4,i)=0.5;
        %Set very low values to J-coupling to avoid it in fitting
        ub(5,i)=MinConstraintRange;
        lb(5,i)=0;
        s0(5,i)=0;        
    elseif strncmp('met', submets{i}, 3)
        %Set boundaries for intensities
        s0(1,i)=lbI + (ubI-lbI)*rand(1);
        ub(1,i)=ubI;
        lb(1,i)=lbI;        
        %Set boundaries for chemical shift
        if constraint_factor(2,i)>0
            %Increments of 0.001
            lbppm=-0.001*constraint_factor(2,i);
            ubppm=0.001*constraint_factor(2,i);
        else
            lbppm=-MinConstraintRange/2;
            ubppm=MinConstraintRange/2;
        end           
        s0(2,i)=fields.pos;        
        ub(2,i)=ubppm+fields.pos;
        lb(2,i)=lbppm+fields.pos;
        %Set boundaries for signal width
        if constraint_factor(3,i)>0
            %Increments of 0.05 Hz (normal ranges 2 Hz)
            lfwhh=(-0.05*constraint_factor(3,i)*1e6)/freq;
            ufwhh=(0.05*constraint_factor(3,i)*1e6)/freq;    
        else
            lfwhh=-MinConstraintRange/2;
            ufwhh=MinConstraintRange/2;
        end
        s0(3,i)=(fields.fwhh*1e6)/freq;
        ub(3,i)=s0(3,i)+ufwhh;
        lb(3,i)=max(s0(3,i)+lfwhh, 0);              
        %Set boundaries for gaussian percentage         
        ugaussFactor=0.1*constraint_factor(4,i)+MinConstraintRange/2;
        lgaussFactor=-0.1*constraint_factor(4,i)-MinConstraintRange/2;
        s0(4,i)=fields.gauss;
        ub(4,i)=min(ugaussFactor+fields.gauss, 1);
        lb(4,i)=max(lgaussFactor+fields.gauss, 0);
        %Set boundaries for J-coupling
        if constraint_factor(5,i)>0,
            %Increments of 0.1 Hz (normal ranges 7 Hz)
            ljcop=(-0.1*constraint_factor(5,i)*1e6)/freq;
            ujcop=(0.1*constraint_factor(5,i)*1e6)/freq;    
        else
            ljcop=0;
            ujcop=MinConstraintRange;
        end
        s0(5,i)=(fields.jcop*1e6)/freq;
        ub(5,i)=s0(5,i)+ujcop;
        lb(5,i)=s0(5,i)+ljcop;                              
    elseif strncmp('spec', submets{i}, 4)
        ub(1:5,i)=MinConstraintRange;
        lb(1:5,i)=0;
        %Set boundaries for intensities
        ub(1,i)=ubI;
        lb(1,i)=lbI;
        s0(:,i)=lb(:,i) +(ub(:,i)-lb(:,i))*rand(1);                
        %If position in fitting for specs doesn´t change is because the
        %change in that variable is so small that it doesn´t move the spectra and then 
        %the optimizer considers that position doesn´t have influence
        %in the fit and stop changing it. One possible solution is
        %resampling, moving and downsampling again        
        n=length(axisscale);
        res=axisscale(1)-axisscale(2);
        ntimes=ceil(max(res/(0.005*1e6/freq), 1)); 
        if ntimes>1,
            xaxisAux=linspace(max(fields.xaxis), min(fields.xaxis), n*ntimes);
            fields.data=interp1(fields.xaxis, fields.data, xaxisAux);
        end
        %First determine the center of the standard
        specCenter=fields.stdPos;
        %Then calculate the different between the center and the initial position
        diffPos=fields.pos-specCenter;
        %Move limits of the standard and create the new axis
        fields.hlim=fields.hlim+diffPos;
        fields.llim=fields.llim+diffPos;
        fields.xaxis=fields.hlim:-(fields.hlim-fields.llim)/(length(fields.data)-1):fields.llim;
        [~, fields.center]=min(abs(fields.xaxis-fields.pos));
        if constraint_factor(2,i)>0
            %Increments of 0.001
            lbppm=-0.001*constraint_factor(2,i);
            ubppm=0.001*constraint_factor(2,i);
        else
            lbppm=-MinConstraintRange/2;
            ubppm=MinConstraintRange/2;
        end           
        s0(2,i)=fields.pos;        
        ub(2,i)=ubppm+fields.pos;
        lb(2,i)=lbppm+fields.pos;
        
        %bRuben
        %Set boundaries for signal width
        if constraint_factor(3,i)~=0
            %Increments of 0.05 Hz (normal ranges 2 Hz)
            lfwhh=(-0.05*constraint_factor(3,i)*1e6)/freq;
            ufwhh=(0.05*constraint_factor(3,i)*1e6)/freq;    
        else
            lfwhh=(-MinConstraintRange/2*1e6)/freq;
            ufwhh=(MinConstraintRange/2*1e6)/freq;
        end
        s0(3,i)=(fields.fwhh*1e6)/freq;
        ub(3,i)=s0(3,i)+ufwhh;
        lb(3,i)=s0(3,i)+lfwhh;                    
        
        %Set boundaries for gaussian percentage         
        ugaussFactor=0.1*constraint_factor(4,i)+MinConstraintRange/2;
        lgaussFactor=-0.1*constraint_factor(4,i)-MinConstraintRange/2;
        s0(4,i)=0;
        ub(4,i)=min(ugaussFactor, 1);
        lb(4,i)=max(lgaussFactor, 0);
        s0(5,i)=0;
        mets.(submets{i})=fields;
    end
end

%Set background functions
if bitand(4, BGtype),
    %Include background functions based on series of cosines
    order=round((max(axisscale)-min(axisscale))/0.015);
    order=min(12, order);
    %Increment counter of the number of variables to fit
    nbnVariables=nbnVariables+order*5;
    %Define cosine functions as additional functions in mets structure
    for i=1:order,
        mets.(['cos' num2str(i)])=[];
    end
    %Define boundaries (necessarily of length 5 to fit mets constraints)
    ub_bg=[1; 2*pi; MinConstraintRange; MinConstraintRange; MinConstraintRange]; %intensity and phase
    lb_bg=[0;    0; 0; 0; 0];
    s0_bg=lb_bg +(ub_bg-lb_bg)*rand(1);
    s0_bg(3:5)=0;
    ub_bg=repmat(ub_bg, 1, order);
    lb_bg=repmat(lb_bg, 1, order);
    s0_bg=repmat(s0_bg, 1, order);
    ub=cat(2, ub, ub_bg);
    lb=cat(2, lb, lb_bg);
    s0=cat(2, s0, s0_bg); 
end
if bitand(2, BGtype),
    minXLim=min(axisscale);
    nbnFunctions=round((max(axisscale)-minXLim)/0.025);
    %Define baseline functions as additional functions in mets structure
    ubAux=zeros(5, nbnFunctions);
    lbAux=zeros(5, nbnFunctions);
    s0Aux=zeros(5, nbnFunctions);
    for i=1:nbnFunctions,
        mets.(['bckg' num2str(i)])=[];
        ubAux(1,i)=0.5;
        lbAux(1,i)=0;  
        s0Aux(1,i)=lbAux(1,i) + (ubAux(1,i)-lbAux(1,i))*rand(1);
        lbAux(2,i)=minXLim+0.04*(i-1);
        ubAux(2,i)=minXLim+0.04*i;        
        s0Aux(2,i)=lbAux(2,i) + (ubAux(2,i)-lbAux(2,i))*rand(1);
        ubAux(3,i)=ufwhhbckg;
        lbAux(3,i)=lfwhhbckg;
        s0Aux(3,i)=lbAux(3,i) + (ubAux(3,i)-lbAux(3,i))*rand(1);
        ubAux(4,i)=ugaussFactorbckg;
        lbAux(4,i)=lgaussFactorbckg;        
        s0Aux(4,i)=0.5;
        %Set very low values to J-coupling to void it in fitting
        ubAux(5,i)=MinConstraintRange;
        lbAux(5,i)=0;
        s0Aux(5,i)=0;        
    end 
    ub=cat(2, ub, ubAux);
    lb=cat(2, lb, lbAux);
    s0=cat(2, s0, s0Aux);     
end
if bitand(1, BGtype),
    mets.poly=[];
    ub_bg=[1; 1; MinConstraintRange; MinConstraintRange; MinConstraintRange]; %offset and slope
    lb_bg=[-1; -1; 0; 0; 0];       
    s0_bg=lb_bg +(ub_bg-lb_bg)*rand(1);
    s0_bg(3:5)=0;
    ub=cat(2, ub, ub_bg);
    lb=cat(2, lb, lb_bg);
    s0=cat(2, s0, s0_bg);         
end

% Setup optimization options
options = optimset('lsqcurvefit'); 
% options = optimset('Algorithm', 'levenberg-marquardt');
options = optimset(options,'DiffMinChange', DiffMinChange);
options = optimset(options,'DiffMaxChange', 1);
% options = optimset(options,'DiffMaxChange', max(sum(ub)));
% options = optimset(options,'DiffMinChange', 1*10^(-1*(7+scale_factor)));
options = optimset(options,'TolFun' ,1*10^(-1*(5+scale_factor))); %Min change allowed in xdata
options = optimset(options,'TolX' ,1*10^(-1*(5+scale_factor))); %Min change allowed in x
options = optimset(options,'MaxIter', 400); %Maximum number of trials with different x0 (reached if stopping creiteria is not fulfill)
options = optimset(options,'MaxFunEvals', nbnVariables*100); %Maximum number of evaluations of TolFun in one iteration (reached if stopping creiteria is not fulfill)
options = optimset(options,'LargeScale' ,'on');
options = optimset(options,'Display' ,'final');

options.maxError = maxError;
options.tries = tries;

Error_best=100;
indexTries=0;

while Error_best>options.maxError && indexTries < options.tries,
    %launch the curve fitting function
    [s,resnorm,~,flag]=lsqcurvefit(@(x, xdata) expModel(x, xdata, mets), s0, axisscale, X, lb, ub, options);    
    
    Error=sqrt(resnorm/length(X))*100/((max(X)-min(X)));
    
    if Error<Error_best,
        Error_best=Error;
        flag_best=flag;
        s_best=s;
    end
    
    if indexTries>=options.tries/2,
        zeroElements=(s0==0);
        for j=1:size(s0,1),
            s0(j,:)=lb(j,:) + (ub(j,:)-lb(j,:)).*rand(1,size(s0,2));
        end
        s0(zeroElements)=0;
    else
        s0(1,:)=lb(1,:) + (ub(1,:)-lb(1,:)).*rand(1,size(s0,2));
    end
    indexTries=indexTries+1;
end

output_args.flag=flag_best;
output_args.F=expModel(s_best, axisscale, mets);
submets=fieldnames(mets);
j=1;

pureSpectra=zeros(length(submets), length(axisscale));

for i=1:length(submets),
    mets_aux=[];
    mets_aux.(submets{i})=mets.(submets{i});
    pureSpectra(i,:,:)=expModel(s_best(:,i), axisscale, mets_aux);
    if strncmp('met', submets{i}, 3) || strncmp('spec', submets{i}, 4),
        s_best([3 5],i)=s_best([3 5],i)*freq/1e6;
        output_args.s(:,j)=s_best(:,i);
        j=j+1;
    end
end
output_args.error=Error_best;

output_args.s(1,:)=output_args.s(1,:)*scale_factor_intensities;
X=X*scale_factor_intensities;
output_args.F=output_args.F*scale_factor_intensities;
pureSpectra=pureSpectra*scale_factor_intensities;

output_args.h=figure;
set(output_args.h, 'visible', 'off');

f(1)=plot(axisscale, X, 'r');
hold on;        
f(2)=plot(axisscale, output_args.F, 'g');
f_back=zeros(1,length(axisscale));
j=1;
for i=1:size(pureSpectra,1),
    if ~strncmp('met', submets{i}, 3) && ~strncmp('spec', submets{i}, 4)
        f_back=f_back+pureSpectra(i,:);
    else
        output_args.pureSpectra(j,:)=pureSpectra(i,:);
        f(2+j)=plot(axisscale, pureSpectra(i,:), 'b'); 
        j=j+1;
    end
end
output_args.baselineSpectra=f_back;
f(end+1)=plot(axisscale, f_back, 'k');
title(['Error: ' num2str(Error_best)]);
set(gca, 'xdir', 'reverse');

if plotFinal,
    set(output_args.h, 'visible', 'on');
end

end

function F = expModel(s, xaxis, mets)
    center=s(2,:);
    intensities=s(1,:);
        
    submets=fieldnames(mets);
    F=0;
    for i=1:length(submets),
        if strncmp('bckg', submets{i}, 4) || strncmp('met', submets{i}, 3)
            fields=getfield(mets, submets{i});
            %First obtain the ratios of intensities depending on the
            %multiplicity
            if strncmp('bckg', submets{i}, 4),
                subPascal=pascal(1);
            else
                subPascal=pascal(fields.mult);
                if size(subPascal, 1)>1,
                    subPascal=subPascal(size(subPascal,1):size(subPascal,1)-1:end-1); 
                end
            end
            %Calcultate the separation of the first signal from the center
            if rem(length(subPascal), 2)==0,
                offset=0.5*s(5,i)+((length(subPascal)/2)-1)*s(5,i);
            else
                offset=floor(length(subPascal)/2)*s(5,i);            
            end
            offset2=offset+center(1,i);
            %generate the j-coupling profile 
            mulPeak=zeros(1, length(xaxis));
            for j=1:length(subPascal),
                individualCenter=offset2-(j-1)*s(5,i);
                peakPars.int=(intensities(1,i)/max(subPascal))*subPascal(j);
                peakPars.pos=individualCenter;
                peakPars.hwhh=s(3,i)/2;
                peakPars.gauss=s(4,i);
                mulPeak = mulPeak + voigtProfile(peakPars, xaxis);
            end
            F=F+mulPeak;
        elseif strncmp('spec', submets{i}, 4),
            %If position in fitting for specs doesn´t change is because the
            %change in that variable is so small that it doesn´t move the spectra and then 
            %the optimizer considers that position doesn´t have influence
            %in the fit and stop changing it. One possible solution is
            %resampling, moving and downsampling again
            spec=mets.(submets{i}).data*s(1,i);
            scaleCorrection=max(max(spec));

            [~, newPos]=min(abs(mets.(submets{i}).xaxis-s(2,i)));
            shift=newPos-mets.(submets{i}).center;
            if shift>0,
                xshift=mets.(submets{i}).xaxis(1)-mets.(submets{i}).xaxis(shift+1);
                addX=linspace(mets.(submets{i}).xaxis(end), mets.(submets{i}).xaxis(end)-xshift, shift+1);
                mets.(submets{i}).xaxis=[mets.(submets{i}).xaxis(shift+1:end) addX(2:end)];
            elseif shift<0
                shift=abs(shift);
                xshift=mets.(submets{i}).xaxis(end-shift)-mets.(submets{i}).xaxis(end);                   
                addX=linspace(mets.(submets{i}).xaxis(1)+xshift, mets.(submets{i}).xaxis(1), shift+1);
                mets.(submets{i}).xaxis=[addX(1:end-1) mets.(submets{i}).xaxis(1:end-shift)];                    
            end
            spec=interp1(mets.(submets{i}).xaxis, spec, xaxis);
            spec(isnan(spec))=0;

            if s(3,i)~=0,
                isnegative=s(3,i)~=abs(s(3,i));
                s(3,i)=abs(s(3,i));
                %Give two passes to avoid asymmetry
                peakPars.int=1;
                peakPars.pos=mean(xaxis);
                peakPars.hwhh=s(3,i)/4;
                peakPars.gauss=s(4,i);
                shapeFactor=voigtProfile(peakPars, xaxis);
                shapeFactor=shapeFactor./sum(shapeFactor);
                if ~isnegative,
                    for pass=1:2,
                        convSpec=conv(spec, shapeFactor);
                        spec=convSpec(ceil(length(shapeFactor)/2):end-floor(length(shapeFactor)/2));
                        spec=fliplr(spec);
                    end
                else
                    %Any ideas?
                end
                spec=spec.*(scaleCorrection/max(max(spec)));
            end
            F=F+spec;
        else
            if strncmp('cos', submets{i}, 3),
                trigoOrder=0.5*str2double(submets{i}(4));
                %Last factor in the equation avoids negative
                cosineBG=s(1,i)*cos(trigoOrder*pi.*(1:length(xaxis))/length(xaxis) + s(2,i))+s(1,i);
                F=F+cosineBG;
            end
            if strncmp('poly', submets{i}, 4),
                coeffs=[s(1,i) s(2,i)];
                f=polyval(coeffs, (0:(length(xaxis)-1))/(length(xaxis)-1));
                F=F+f;
            end
        end
    end
end