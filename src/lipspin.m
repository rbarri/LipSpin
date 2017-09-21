function lipspin()
%   LipSpin() - 1H-NMR lipid profiling GUI
%   Please consult the supplied manual for general information about the
%   program
%
%   Example: Start the graphical user interface by typing "Lipspin" in the
%            command window provided that the application folder is
%            included in the general search path

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

%Initialise
NMRdata.version='1.0';
NMRdata.phaseRegions=[];
NMRdata.bcRegions=[];

%if exists then read config file, otherwise set default values 
if exist('options.nmrcfg', 'file') == 2,
    f=fopen('options.nmrcfg');
    cfgFile = fread(f, '*char')';
    fclose(f);
    if ~strcmp(cfgFile,'')
        line=1;
        posCR=unique([0 strfind(cfgFile, sprintf('\n')) length(cfgFile)]);
        while (line<length(posCR))
            currentLine=cfgFile(posCR(line)+1:posCR(line+1));
            aux=strfind(currentLine, '=');
            switch currentLine(3:aux(1)-1)          
                case 'MAXERROR'
                    NMRdata.fitPars.maxError=str2double(currentLine(aux+1:end-1));
                case 'MAXITERS'
                    NMRdata.fitPars.maxIters=str2double(currentLine(aux+1:end-1));
                case 'NORMACQPARS'
                    NMRdata.fitPars.normAcqPars=str2double(currentLine(aux+1:end-1));
                case 'PATTERNS_FOLDER'
                    NMRdata.paths.patterns=(currentLine(aux+1:end-1));
                case 'STDS_FOLDER'
                    NMRdata.paths.stds=(currentLine(aux+1:end-1));
            end
            line=line+1;
        end
    end
else
    %Stop condition for fitting. Maximum error allowed in %
    NMRdata.fitPars.maxError=1;
    %Maximum number of fitting attempts without reaching stop condition
    NMRdata.fitPars.maxIters=3;
    %Normalise spectra by RG, NS and P1
    NMRdata.fitPars.normAcqPars=1;
    NMRdata.paths.patterns=pwd;
    NMRdata.paths.stds=pwd;    
end

%%  Main Window Setup
%Figure and axes
set(0,'units','pixels');  
NMRdata.screenRes = get(0,'screensize');
hMainFigure = figure('Units','pixels', 'MenuBar','none', 'Name', ['LipSpin v.' NMRdata.version], 'NumberTitle','Off', 'DeleteFcn',@hMainFigure_DeleteFcn, 'tag', 'MainFigure', 'Toolbar','Figure', 'OuterPosition',[0.0 0.0 NMRdata.screenRes(3) NMRdata.screenRes(4) ], 'Visible','off');
personaliseToolbar(hMainFigure);
NMRdata.hMainAxes = axes('Parent', hMainFigure, 'Units','normalized', 'Position',[0.02 0.5 0.85 0.45]);
aux=get(hMainFigure, 'Color');

%UIcontrols
%axis labels
NMRdata.lblMainXaxis=uicontrol(hMainFigure, 'Style','text', 'Units','Normalized', 'backgroundColor', aux, 'HorizontalAlignment','left', 'Position',[0.4 0.42 0.1 0.05], 'String', 'chemical shift (ppm)', 'ButtonDownFcn', @lblMainXaxis_Callback ); 
%included & visible samples
pnlSamples=uipanel(hMainFigure,'Units','Normalized', 'Position',[0.88 0.5 0.11 0.45]);
NMRdata.uiControls.btnSelAll=uicontrol(pnlSamples, 'Style','PushButton', 'Units','Normalized', 'Position',[0.05 0.02 0.5 0.06], 'String','Select All', 'CallBack', {@btnSelAll_Callback}, 'enable', 'off');
NMRdata.uiControls.tblIncludedSamples = uitable(pnlSamples, 'Units', 'Normalized', 'RowName', [], 'Position', [0.05 0.1 0.90 0.45], 'Data', [], 'ColumnName', {'Include', 'Sample'}, 'ColumnEditable', [true false], 'CellEditCallback', @tblIncludedSamples_CellEditCallback);
set(NMRdata.uiControls.tblIncludedSamples, 'ColumnFormat', {'logical', 'numeric'});
lblIncluded=uicontrol(pnlSamples, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.57 0.85 0.04], 'String', 'Included sample(s)'); 
NMRdata.uiControls.lbVisible=uicontrol(pnlSamples, 'Style','listbox', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.05 0.63 0.90 0.3], 'String', [], 'max',2, 'min',0, 'CallBack', {@lbVisible_Callback}, 'enable', 'off');
lblVisible=uicontrol(pnlSamples, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.94 0.85 0.04], 'String', 'Visible sample(s)'); 
%ROI list
pnlRegionList=uipanel(hMainFigure, 'FontWeight','bold', 'TitlePosition','lefttop', 'title', 'Signal patterns','Units','Normalized', 'Position',[0.02 0.02 0.18 0.4]);
btnLoadPattern=uicontrol(pnlRegionList, 'Style','PushButton', 'Units','Normalized', 'Position',[0.1 0.9 0.35 0.08 ], 'String','Load pattern', 'CallBack', {@btnLoadPattern_Callback});        
NMRdata.uiControls.btnSavePattern=uicontrol(pnlRegionList, 'Style','PushButton', 'Units','Normalized', 'Position',[0.55 0.9 0.35 0.08 ], 'String','Save pattern', 'CallBack', {@btnSavePattern_Callback}, 'Enable', 'off');        
btnCreatePattern=uicontrol(pnlRegionList, 'Style','PushButton', 'Units','Normalized', 'Position',[0.65 0.1 0.25 0.08 ], 'String','Create', 'CallBack', {@btnCreatePattern_Callback});        
NMRdata.uiControls.txtCreateRegion=uicontrol(pnlRegionList, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.1 0.1 0.5 0.08], 'String','');        
NMRdata.uiControls.lbPatterns=uicontrol(pnlRegionList, 'Style','listbox', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.1 0.25 0.5 0.5], 'String', [], 'Callback', @lbPatterns_Callback);        
NMRdata.uiControls.btnDelPattern=uicontrol(pnlRegionList, 'Style','PushButton', 'Units','Normalized', 'Position',[0.65 0.25 0.25 0.08 ], 'String','Remove', 'CallBack', {@btnDelPattern_Callback}, 'Enable', 'off');        
%ROI properties
NMRdata.pnlRegionSettings=uipanel(hMainFigure, 'FontWeight','bold', 'TitlePosition','lefttop', 'title', 'Pattern Settings','Units','Normalized', 'Position',[0.21 0.02 0.2 0.4]);
lblMaxPPM=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.1 0.9 0.2 0.05], 'String', 'Max ppm: ');        
NMRdata.uiControls.txtMaxPPM=uicontrol(NMRdata.pnlRegionSettings, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.9 0.15 0.05], 'String','', 'Callback', @lblMaxPPM_Callback);
lblMinPPM=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.55 0.9 0.2 0.05], 'String', 'Min ppm: ');         
NMRdata.uiControls.txtMinPPM=uicontrol(NMRdata.pnlRegionSettings, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.75 0.9 0.15 0.05], 'String','', 'Callback', @lblMinPPM_Callback);        
lblFitMode=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.1 0.8 0.35 0.05], 'String','Select Mode: ' ); 
NMRdata.uiControls.popFitMode=uicontrol(NMRdata.pnlRegionSettings, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.55 0.8 0.35 0.05], 'String', {'Lineshape fitting', 'Integration'}, 'Callback',@popFitMode_Callback);                 
lblFitBaseline=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.1 0.7 0.5 0.05], 'String','Baseline signals' ); 
lblFitBaselinePol=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.1 0.63 0.25 0.05], 'String','1st Poly: ' ); 
NMRdata.uiControls.cbFitBaselinePol=uicontrol(NMRdata.pnlRegionSettings, 'Style','checkbox', 'Units','Normalized', 'Position',[0.33 0.63 0.045 0.05], 'Callback', @cbFitBaseline_Callback);                
lblFitBaselineGauss=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.44 0.63 0.18 0.05], 'String','Gaussians: ' ); 
NMRdata.uiControls.cbFitBaselineGauss=uicontrol(NMRdata.pnlRegionSettings, 'Style','checkbox', 'Units','Normalized', 'Position',[0.62 0.63 0.045 0.05], 'Callback', @cbFitBaseline_Callback);  
lblFitBaselineCos=uicontrol(NMRdata.pnlRegionSettings, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.73 0.63 0.12 0.05], 'String','Cosine: ' ); 
NMRdata.uiControls.cbFitBaselineCos=uicontrol(NMRdata.pnlRegionSettings, 'Style','checkbox', 'Units','Normalized', 'Position',[0.85 0.63 0.045 0.05], 'Callback', @cbFitBaseline_Callback);  
NMRdata.uiControls.lbMetList=uicontrol(NMRdata.pnlRegionSettings, 'Style','listbox', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.1 0.2 0.5 0.4], 'String', [], 'Callback', @lbMetList_Callback);        
NMRdata.uiControls.btnDelMet=uicontrol(NMRdata.pnlRegionSettings, 'Style','PushButton', 'Units','Normalized', 'Position',[0.65 0.2 0.25 0.08 ], 'String','Delete', 'CallBack', {@btnDelMet_Callback});        
NMRdata.uiControls.btnNewMet=uicontrol(NMRdata.pnlRegionSettings, 'Style','PushButton', 'Units','Normalized', 'Position',[0.65 0.1 0.25 0.08 ], 'String','New', 'CallBack', {@btnNewMet_Callback});             
NMRdata.uiControls.txtNewMet=uicontrol(NMRdata.pnlRegionSettings, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.1 0.1 0.5 0.08], 'String','');        
h=get(NMRdata.pnlRegionSettings, 'children');
EnableDisableUIcontrols(h, 'off');
%Additional fitting parameters
lblCorrectWidth=uicontrol(hMainFigure, 'Style','text', 'Units','Normalized', 'fontweight', 'bold','backgroundColor', aux, 'HorizontalAlignment','left', 'Position', [0.43 0.38 0.1 0.02], 'String', 'Correct signal width by:');        
NMRdata.uiControls.popCorrectWidth=uicontrol(hMainFigure, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.43 0.35 0.1 0.02], 'String', {'None'}, 'Callback',@popCorrectWidth_Callback);                 
NMRdata.uiControls.lblRefSample=uicontrol(hMainFigure, 'Style','text', 'Units','Normalized', 'fontweight', 'bold','backgroundColor', aux, 'HorizontalAlignment','left', 'visible', 'off', 'Position', [0.43 0.31 0.1 0.02], 'String', 'Reference sample:');        
NMRdata.uiControls.popRefSample=uicontrol(hMainFigure, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.43 0.28 0.1 0.02], 'visible', 'off', 'String', '', 'Callback', @popRefSample_Callback);                 
lblFigPath=uicontrol(hMainFigure, 'Style','text', 'Units','Normalized', 'fontweight', 'bold','backgroundColor', aux, 'HorizontalAlignment','left', 'Position', [0.43 0.15 0.1 0.02], 'String', 'Figures path');        
NMRdata.uiControls.txtFigPath=uicontrol(hMainFigure, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.43 0.12 0.1 0.03], 'String','');
cmPath = uicontextmenu('Parent', hMainFigure);
mnPath = uimenu(cmPath,'Label','paste', 'CallBack',{@mnPath_Callback, NMRdata.uiControls.txtFigPath});
set(NMRdata.uiControls.txtFigPath, 'UIContextMenu', cmPath);
btnFit=uicontrol(hMainFigure, 'Style','PushButton', 'Units','Normalized', 'fontsize', 14, 'Position',[0.43 0.02 0.1 0.07 ], 'String','Run', 'CallBack', {@btnFit_Callback});             
%Results panel
pnlFitResults=uipanel(hMainFigure, 'FontWeight','bold', 'TitlePosition','lefttop', 'title', 'Results','Units','Normalized', 'Position',[0.55, 0.02, 0.44, 0.40]);
btnResultLastPar=uicontrol(pnlFitResults, 'Style','PushButton', 'Units','Normalized', 'Fontsize', 14, 'Position',[0.025 0.9 0.05 0.08 ], 'String','<', 'CallBack', {@btnParResult_Callback, -1});             
NMRdata.uiControls.lblParResult=uicontrol(pnlFitResults, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','center', 'Position', [0.075 0.89 0.15 0.1 ], 'String', ' ');        
btnResultNextPar=uicontrol(pnlFitResults, 'Style','PushButton', 'Units','Normalized', 'Fontsize', 14, 'Position',[0.225 0.9 0.05 0.08 ], 'String','>', 'CallBack', {@btnParResult_Callback,  1});             
btnResultExport=uicontrol(pnlFitResults, 'Style','PushButton', 'Units','Normalized', 'Position',[0.8 0.9 0.175 0.08 ], 'String','Export table', 'CallBack', {@btnResultExport_Callback});             
btnResultDelCol=uicontrol(pnlFitResults, 'Style','PushButton', 'Units','Normalized', 'Position',[0.6 0.9 0.19 0.08 ], 'String','Remove selected column', 'CallBack', {@btnResultDelCol_Callback});             
for i=1:6
    NMRdata.uiControls.tblresults(i) = uitable(pnlFitResults, 'Units', 'Normalized', 'RowName', 'sample', 'visible', 'off', 'Position', [0.025, 0.1, 0.95, 0.75], 'RowName',[], 'Columnname',[], 'Data', [], 'CellSelectionCallback', {@tblresults_CellSelectionCallback});
end
set(NMRdata.uiControls.tblresults(1), 'visible', 'on');

NMRdata.fitPars.correctShape.sampleRef=1;
NMRdata.fitPars.correctShape.widthRef='None';

% Menu bar
%File menu
hMenuFile=uimenu(hMainFigure, 'Label','File');
% uimenu(hMenuFile,'Label','Open', 'Enable','On', 'Callback',@Open_data);
hImportFIDMenu=uimenu(hMenuFile,'Label','Import FID', 'Enable','On');
uimenu(hImportFIDMenu,'Label','Bruker', 'Enable','On', 'tag', 'Bruker FID', 'Callback',@ImportBruker_Callback);
hImportSpecMenu=uimenu(hMenuFile,'Label','Import spectra', 'Enable','On');
uimenu(hImportSpecMenu,'Label','Bruker', 'Enable','On', 'tag', 'Bruker spectra', 'Callback',@ImportBruker_Callback);
%Pre-processing menu
hMenuProc=uimenu( hMainFigure, 'Label','Pre-processing');
hMenuAligner=uimenu(hMenuProc,'Label','Align', 'Enable','On', 'Callback', @Align_Callback);
hMenuPhase=uimenu(hMenuProc,'Label','Autophase', 'Enable','On', 'Callback', @Phase_Callback);
hMenuBaseline=uimenu(hMenuProc,'Label','Baseline', 'Enable','On', 'Callback', @Baseline_Callback);
hMenuRefShift=uimenu(hMenuProc,'Label','Reference Chemical Shift', 'Enable','On', 'Callback', @RefShift_Callback);
hMenuRefDeconv=uimenu(hMenuProc,'Label','Reference Deconvolution', 'Enable','On', 'Callback', @RefDeconv_Callback);
%Standard profiles menu
hMenuStd=uimenu( hMainFigure, 'Label','Standard spectra'); 
hMenuSaveStandard=uimenu(hMenuStd,'Label','Save standard', 'Enable','On', 'Callback', @SaveStandard_Callback);
hMenuPlotStandard=uimenu(hMenuStd,'Label','Load session standards', 'Enable','On', 'Callback', @LoadStandards_Callback);
%General options menu
hMenuOptions=uimenu( hMainFigure, 'Label','Options');
hMenuGeneralOptions=uimenu(hMenuOptions,'Label','General Options', 'Enable','On', 'Callback', @GeneralOptions_Callback);

%Switch main screen visible
set(hMainFigure,'Visible','on')

%%  Exiting functions 
    function hMainFigure_DeleteFcn(~,~)  
        QuitLipspin();
    end

    function hImportOptions_DeleteFcn(~,~)
        NMRdata=rmfield(NMRdata,'hImportOptions');
        if isfield(NMRdata,'importData'),
            NMRdata=rmfield(NMRdata,'importData');
        end
    end

    function hPhaseOptions_DeleteFcn(~,~)
        NMRdata=rmfield(NMRdata,'hPhaseOptions'); 
        if isfield(NMRdata,'phaseData'),
            NMRdata=rmfield(NMRdata,'phaseData');
        end
    end

    function hBaselineOptions_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hBaselineOptions');
        if isfield(NMRdata,'baselineData'),
            NMRdata=rmfield(NMRdata,'baselineData');
        end
    end

    function hAlignOptions_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hAlignOptions'); 
        if isfield(NMRdata,'alignData'),
            NMRdata=rmfield(NMRdata,'alignData');
        end
    end

    function hRefShiftOptions_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hRefShiftOptions');   
        if isfield(NMRdata, 'refShiftData')
            NMRdata=rmfield(NMRdata,'refShiftData');   
        end          
    end

    function hRefDeconvOptions_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hRefDeconvOptions');  
        if isfield(NMRdata, 'refDeconvData')
            NMRdata=rmfield(NMRdata,'refDeconvData');
        end  
    end

    function hEditGeneralOptions_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hEditGeneralOptions');
        if isfield(NMRdata, 'generalOptionsData')
            NMRdata=rmfield(NMRdata, 'generalOptionsData');
        end
    end

    function hLoadStandards_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hLoadStandards');   
        if isfield(NMRdata, 'loadStandardsData')
            NMRdata=rmfield(NMRdata, 'loadStandardsData');
        end
    end

    function hEditMet_DeleteFcn (~,~)
        NMRdata=rmfield(NMRdata,'hEditMet');   
        if isfield(NMRdata, 'editMetData')
            NMRdata=rmfield(NMRdata,'editMetData');   
        end
    end

%% Menu bar callbacks
    function ImportBruker_Callback(hObject, ~)
        if isfield(NMRdata,'hImportOptions')
            close(NMRdata.hImportOptions);
        end
        mode=get(hObject, 'tag');
        
        NMRdata.hImportOptions= figure('Units','pixels', 'Name', 'Import Options', 'NumberTitle', 'off', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hImportOptions_DeleteFcn, 'Visible','off');

        pnlOptionsPanel=uipanel('Parent',NMRdata.hImportOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);

        lblPath=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.85 0.3 0.05], 'String','Path files: ' ); 
        NMRdata.importData.txtPath=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.85 0.2 0.05], 'String','');
        cmPath = uicontextmenu('Parent', NMRdata.hImportOptions);
        mnPath = uimenu(cmPath,'Label','paste', 'CallBack',{@mnPath_Callback, NMRdata.importData.txtPath});
        set(NMRdata.importData.txtPath, 'UIContextMenu', cmPath);
        btnPath=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.81 0.85 0.1 0.05 ], 'String','Search', 'CallBack', {@searchDirPath_Callback});         
        lblFilter=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.75 0.3 0.05], 'String','File name filter: ' ); 
        NMRdata.importData.txtFilter=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.75 0.2 0.05], 'String','');                
        lblExp_no=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.65 0.3 0.05], 'String','Experiment number: ' ); 
        NMRdata.importData.txtExp_no=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.65 0.2 0.05], 'String','');
        switch mode
            case 'Bruker spectra'
                NMRdata.importData.typeFile=1;
                lblProc_no=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.55 0.3 0.05], 'String','Processing number: ' ); 
                NMRdata.importData.txtProc_no=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.55 0.2 0.05], 'String','');
            case 'Bruker FID'
                NMRdata.importData.typeFile=2;
                NMRdata.importData.WDWtype='non';
                lblZerofilling=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.5 0.3 0.05], 'String','Zero filling factor in FT: ' ); 
                lblZerofillingF2=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.77 0.55 0.05 0.05], 'String','F2' );                 
                NMRdata.importData.txtZerofillingF2=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.72 0.5 0.08 0.05], 'String','0');     
                lblFirstPointby05=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.2 0.4 0.3 0.05], 'String','Multiply first point by 0.5: ' ); 
                NMRdata.importData.cbFirstPointby05=uicontrol('Parent',pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized', 'Position',[0.774 0.4 0.04 0.04]);                
                pnlOptionsWDW=uipanel('Parent',pnlOptionsPanel, 'TitlePosition','lefttop', 'Units','Normalized', 'Position',[0.2 0.1 0.5 0.25], 'title', 'Apodization options');
                lblWDWF2=uicontrol('Parent',pnlOptionsWDW, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.65 0.3 0.2], 'String','Type: ' ); 
                popWDWF2=uicontrol('Parent',pnlOptionsWDW, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.55 0.65 0.3 0.25], 'String', {'none', 'lorentzian', 'gaussian'}, 'Callback', @popWDW_Callback);                
                lblWDWF2lb=uicontrol('Parent',pnlOptionsWDW, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.15 0.2 0.2], 'String','LB: ' ); 
                NMRdata.importData.txtWDWF2lb=uicontrol('Parent',pnlOptionsWDW, 'Style','edit', 'Units','Normalized', 'backgroundColor',' w', 'Position',[0.25 0.15 0.2 0.25], 'String', '', 'enable', 'off'); 
                lblWDWF2gb=uicontrol('Parent',pnlOptionsWDW, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.55 0.15 0.2 0.2], 'String','GB: ' ); 
                NMRdata.importData.txtWDWF2gb=uicontrol('Parent',pnlOptionsWDW, 'Style','edit', 'Units','Normalized', 'backgroundColor',' w', 'Position',[0.75 0.15 0.2 0.25], 'String', '', 'enable', 'off');               
        end
        btnImport=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.725 0.1 0.2 0.1 ], 'String','Import', 'CallBack', {@btnImport_Callback});        
        
        set(NMRdata.hImportOptions,  'Visible', 'on');
    end
    
    function Phase_Callback(~,~)
        if isfield(NMRdata,'hPhaseOptions')
            close(NMRdata.hPhaseOptions)
        end
        NMRdata.hPhaseOptions=figure('Units','pixels', 'Name', 'Edit phase settings', 'NumberTitle', 'off', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hPhaseOptions_DeleteFcn, 'Visible','off');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hPhaseOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);                        
        lblPhaseMode=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.7 0.2 0.1], 'String','Select method: ' ); 
        NMRdata.phaseData.popPhaseMode=uicontrol('Parent',pnlOptionsPanel, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.4 0.7 0.25 0.1], 'String', {'entropy', 'flat baseline'}, 'CallBack', {@popPhaseMode_Callback});                 
        lblPhaseAutoRegions=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.65 0.45 0.05], 'String','Automatic range detection: ' ); 
        NMRdata.phaseData.cbPhaseAutoRegions=uicontrol('Parent', pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized','Position',[0.62 0.65 0.04 0.04], 'tag', 'phaseData', 'Enable', 'off', 'CallBack', {@cbAutoRange_Callback});               
        NMRdata.phaseData.tblRanges=uitable(pnlOptionsPanel, 'Units','Normalized', 'RowName', [], 'Position',[0.15 0.15 0.5 0.4], 'Data',zeros(1, 2), 'ColumnName',{'higher ppm', 'lower ppm'}, 'ColumnEditable', [true true], 'tag', 'phaseData', 'CellSelectionCallback', {@tblRanges_CellSelectionCallback});
        set(NMRdata.phaseData.tblRanges, 'Units', 'pixel');
        aux=get(NMRdata.phaseData.tblRanges, 'Position');
        set(NMRdata.phaseData.tblRanges,'ColumnWidth', {aux(3)/2  aux(3)/2});
        set(NMRdata.phaseData.tblRanges, 'Units', 'Normalized');
        if ~isempty(NMRdata.phaseRegions),
            set(NMRdata.phaseData.tblRanges, 'Data', NMRdata.phaseRegions);
        end
        NMRdata.phaseData.btnAddRange=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.7 0.5 0.2 0.05 ], 'String','Add region', 'tag', 'phaseData', 'CallBack', {@btnAddRange_Callback});        
        NMRdata.phaseData.btnDelRange=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.7 0.4 0.2 0.05 ], 'String','Remove selected', 'tag', 'phaseData', 'CallBack', {@btnDelRange_Callback});
        NMRdata.phaseData.btnPhase=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.7 0.15 0.2 0.1 ], 'String','Phase', 'CallBack', {@btnPhase_Callback});        
        set(NMRdata.hPhaseOptions,  'Visible', 'on');
    end

    function Baseline_Callback(~,~)     
        if isfield(NMRdata,'hBaselineOptions')
            close(NMRdata.hBaselineOptions)
        end
        NMRdata.hBaselineOptions=figure('Units','pixels', 'Name','Baseline correction settings', 'NumberTitle', 'off', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hBaselineOptions_DeleteFcn, 'Visible','off');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hBaselineOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);                        
        lblBaselineMode=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.75 0.2 0.05], 'String','Select Mode: ' ); 
        NMRdata.baselineData.popBaselineMode=uicontrol('Parent',pnlOptionsPanel, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.4 0.75 0.25 0.05], 'String', {'median substraction', 'cubic spline', 'cubic Hermite', 'polynomial'}, 'CallBack', {@popBaselineMode_Callback});                 
        NMRdata.baselineData.lblOrder=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.7 0.75 0.1 0.05], 'Visible', 'off', 'String', 'Order: ' );    
        NMRdata.baselineData.txtOrder=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.8 0.75 0.1 0.05], 'Visible', 'off', 'String','');
        lblBaselineAutoRegions=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.65 0.45 0.05], 'String','Automatic range detection: ' ); 
        NMRdata.baselineData.cbBaselineAutoRanges=uicontrol('Parent', pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized','Position',[0.62 0.65 0.04 0.04], 'tag', 'baselineData', 'CallBack', {@cbAutoRange_Callback});               
        NMRdata.baselineData.tblRanges=uitable(pnlOptionsPanel, 'Units','Normalized', 'Position',[0.15 0.2 0.5 0.4], 'Data',zeros(1, 2), 'ColumnName',{'higher ppm', 'lower ppm'}, 'ColumnEditable', [true true], 'tag', 'baselineData', 'CellSelectionCallback', {@tblRanges_CellSelectionCallback});
        set(NMRdata.baselineData.tblRanges, 'Units', 'pixel');
        aux=get(NMRdata.baselineData.tblRanges, 'Position');
        set(NMRdata.baselineData.tblRanges,'ColumnWidth', {aux(3)/2.31  aux(3)/2.31});
        set(NMRdata.baselineData.tblRanges, 'Units', 'Normalized');
        if ~isempty(NMRdata.bcRegions),
            set(NMRdata.baselineData.tblRanges, 'Data', NMRdata.bcRegions);
        end        
        NMRdata.baselineData.btnAddRange=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.7 0.55 0.2 0.05 ], 'String','Add region', 'tag', 'baselineData', 'CallBack', {@btnAddRange_Callback});        
        NMRdata.baselineData.btnDelRange=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.7 0.45 0.2 0.05 ], 'String','Remove selected', 'tag', 'baselineData', 'CallBack', {@btnDelRange_Callback});
        NMRdata.baselineData.btnBaseline=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.7 0.2 0.2 0.1 ], 'String','Correct baseline', 'CallBack', {@btnBaseline_Callback});        
        set(NMRdata.hBaselineOptions,  'Visible', 'on');
    end

    function Align_Callback(~,~)
        if isfield(NMRdata,'hAlignOptions')
            close(NMRdata.hAlignOptions)
        end
        NMRdata.hAlignOptions=figure('Units','pixels', 'Name', 'Align settings','NumberTitle', 'off', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hAlignOptions_DeleteFcn, 'Visible','off');
%         aux=get(hAlignOptions, 'Color');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hAlignOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);
        lblAlignRegion=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized','FontWeight','bold','HorizontalAlignment','left', 'Position',[0.15 0.7 0.5 0.05], 'String', 'Introduce region to align: ' );    
        lblMaxPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.6 0.2 0.05], 'String', 'Max ppm: ' );         
        NMRdata.alignData.txtMaxAlignRegion=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.6 0.1 0.05], 'String','');
        lblMinPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.55 0.6 0.2 0.05], 'String', 'Min ppm: ' );         
        NMRdata.alignData.txtMinAlignRegion=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.75 0.6 0.1 0.05], 'String','');        
        lblApplyAllSpectrum=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.44 0.45 0.05], 'String','Shift all spectrum [X] / Only region [ ] : ' ); 
        NMRdata.alignData.cbApplyAllSpectrum=uicontrol('Parent', pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized', 'Position',[0.81 0.45 0.04 0.04]);                
        btnAlign=uicontrol('Parent', pnlOptionsPanel, 'Style', 'PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.65 0.2 0.2 0.1 ], 'String','Align', 'CallBack', {@btnAlign_Callback});                
        set(NMRdata.hAlignOptions,  'Visible', 'on');
    end

    function RefShift_Callback(~,~)
        if isfield(NMRdata,'hRefShiftOptions')
            close(NMRdata.hRefShiftOptions)
        end
        NMRdata.hRefShiftOptions=figure('Units','pixels', 'NumberTitle', 'off', 'Name', 'Reference chemical shift', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hRefShiftOptions_DeleteFcn, 'Visible','off');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hRefShiftOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);        
        lblRefShifCenter=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.25 0.65 0.3 0.05], 'String', 'Chemical shift reference: ' );         
        NMRdata.refShiftData.txtRefShifCenter=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.65 0.15 0.05], 'String','');
        lblRefShifTolerance=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.25 0.55 0.3 0.05], 'String', ['Tolerance (' char(177) ' ppm) : ']);
        NMRdata.refShiftData.txtRefShifTolerance=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.6 0.55 0.15 0.05], 'String','');
        btnRefShif=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.55 0.3 0.2 0.1 ], 'String', 'reference', 'CallBack', {@btnRefShif_Callback});        
        
        set(NMRdata.hRefShiftOptions,  'Visible', 'on');
    end

    function RefDeconv_Callback(~,~)
        if isfield(NMRdata,'hRefDeconvOptions')
            close(NMRdata.hRefDeconvOptions)
        end 
        NMRdata.hRefDeconvOptions=figure('Units','pixels', 'NumberTitle', 'off', 'Name', 'Reference deconvolution settings', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hRefDeconvOptions_DeleteFcn, 'Visible','off');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hRefDeconvOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);                
        lblRefDeconvMode=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.75 0.35 0.05], 'String','Select reference signal: ' ); 
        NMRdata.refDeconvData.popRefDeconvMode=uicontrol('Parent',pnlOptionsPanel, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.5 0.75 0.35 0.05], 'String', {'Signal in dataset', 'Synthetic TMS signal'});                 
        lblRefDeconvRanges=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'FontWeight','bold','HorizontalAlignment','left', 'Position',[0.15 0.62 0.5 0.05], 'String', 'Introduce boundaries: ' );    
        lblRefDeconvMaxPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.55 0.2 0.05], 'String', 'Max ppm: ' );         
        NMRdata.refDeconvData.txtRefDeconvMaxPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.55 0.1 0.05], 'String','');
        lblRefDeconvMinPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.55 0.55 0.2 0.05], 'String', 'Min ppm: ' );         
        NMRdata.refDeconvData.txtRefDeconvMinPPM=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.75 0.55 0.1 0.05], 'String',''); 
        lblRefDeconvSpectrum=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.15 0.4 0.35 0.05], 'String', 'Reference spectrum: ');        
        NMRdata.refDeconvData.popRefDeconvSpectrum=uicontrol('Parent', pnlOptionsPanel, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.5 0.4 0.35 0.05], 'String', NMRdata.X.label{1}); 
        btnRefDeconvSpectrum=uicontrol('Parent', pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'Position',[0.5 0.3 0.35 0.08], 'String', 'Best (width and skewness)', 'CallBack', {@btnRefDeconvSpectrum_Callback}); 
        btnRefDeconv=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.60 0.1 0.25 0.1 ], 'String', 'correct lineshape', 'CallBack', {@btnRefDeconv_Callback});        
      
        set(NMRdata.hRefDeconvOptions,  'Visible', 'on');
    end

    function GeneralOptions_Callback(~,~)
        if isfield(NMRdata,'hEditGeneralOptions')
            close(NMRdata.hEditGeneralOptions)
        end
        NMRdata.hEditGeneralOptions=figure('Units','pixels', 'NumberTitle', 'off', 'Name', 'General Options', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hEditGeneralOptions_DeleteFcn, 'Visible','off');
        pnlOptionsPanel=uipanel('Parent',NMRdata.hEditGeneralOptions, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);        
        lblMaxError=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.1 0.8 0.5 0.05], 'Visible', 'on', 'String', 'Max. error (%) stopping criterion: ' );    
        NMRdata.generalOptionsData.txtMaxError=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.7 0.8 0.2 0.05], 'Visible', 'on', 'String', num2str(NMRdata.fitPars.maxError));
        lblmaxIters=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.1 0.7 0.5 0.05], 'Visible', 'on', 'String', 'Max. number of iterations: ' );    
        NMRdata.generalOptionsData.txtmaxIters=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.7 0.7 0.2 0.05], 'Visible', 'on', 'String', num2str(NMRdata.fitPars.maxIters));
        lblPathStd=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.1 0.6 0.5 0.05], 'Visible', 'on', 'String', 'Default path for standard spectra: ' );    
        NMRdata.generalOptionsData.txtPathStd=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.7 0.6 0.2 0.05], 'Visible', 'on', 'String', strrep(NMRdata.paths.stds, '\\', '\'));
        lblPathPatterns=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.1 0.5 0.5 0.05], 'Visible', 'on', 'String', 'Default path for signal patterns: ' );    
        NMRdata.generalOptionsData.txtPathPatterns=uicontrol('Parent',pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.7 0.5 0.2 0.05], 'Visible', 'on', 'String', strrep(NMRdata.paths.patterns, '\\', '\'));
        lblNormACQPars=uicontrol('Parent', pnlOptionsPanel, 'Style','text', 'Units','Normalized','HorizontalAlignment','left', 'Position',[0.1 0.4 0.75 0.05], 'Visible', 'on', 'String', 'Normalise imported spectra by acquisition parameters (NS, RG, P1)?: ' );    
        NMRdata.generalOptionsData.cbNormACQPars=uicontrol('Parent', pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized', 'Position',[0.87 0.4 0.04 0.04], 'value', NMRdata.fitPars.normAcqPars);                
        btnSaveGeneralOptions=uicontrol('Parent',pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.7 0.2 0.2 0.1 ], 'String', 'Save', 'CallBack', {@btnSaveGeneralOptions_Callback});         
        set(NMRdata.hEditGeneralOptions,  'Visible', 'on');
    end

    function SaveStandard_Callback(~,~)
        if length(find(NMRdata.X.include{1}))==1,
            [fileName, filepath] = uiputfile('*.nmrstd', 'Save file', NMRdata.X.label{1}{NMRdata.X.include{1}});            
            str=['##NAME=' NMRdata.X.label{1}{NMRdata.X.include{1}}  '\n'];
            str=[str '##LOWLIMIT=' num2str(NMRdata.X.axisscale{2}(end))  '\n'];
            str=[str '##HIGHLIMIT=' num2str(NMRdata.X.axisscale{2}(1))  '\n'];
            str=[str '##FREQUENCY=' num2str(NMRdata.X.SF(1))  '\n'];
            str=[str '##SPECTRUM=\n'];
            f = fopen([filepath fileName], 'w+');       
            fprintf(f, str);
            fseek(f, 0, 'eof');
            spc=real(NMRdata.X.data(NMRdata.X.include{1}, :));            
            %normalise to unit sum
            fwrite(f, spc./sum(spc), 'double');
            fclose(f);
        else
            warndlg('Only one spectrum has to be included at a time ')
        end
    end

    function LoadStandards_Callback(~,~)
        if isfield(NMRdata,'hLoadStandards')
            close(NMRdata.hLoadStandards)
        end
        NMRdata.hLoadStandards=figure('Units','pixels', 'Name', 'Session standards', 'NumberTitle', 'off', 'MenuBar','none', 'tag', 'StdFigure', 'Toolbar','figure', 'DeleteFcn',@hLoadStandards_DeleteFcn, 'Visible','off');
        personaliseToolbar(NMRdata.hLoadStandards);
        NMRdata.loadStandardsData.hStandardAxes = axes('Parent', NMRdata.hLoadStandards, 'Units','normalized', 'Position',[0.05 0.35 0.9 0.6]);
        NMRdata.loadStandardsData.lbStandardList=uicontrol(NMRdata.hLoadStandards, 'Style','listbox', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.05 0.05 0.65 0.25], 'max',2, 'min',0, 'String', [], 'CallBack', {@lbStandardList_Callback});        
        btnLoadStandard=uicontrol('Parent',NMRdata.hLoadStandards, 'Style','PushButton', 'Units','Normalized', 'Position',[0.75 0.2 0.2 0.08 ], 'String','Load Standard', 'CallBack', {@btnLoadFileStandard_Callback});        
        NMRdata.loadStandardsData.btnRemoveStandard=uicontrol('Parent',NMRdata.hLoadStandards, 'Style','PushButton', 'Units','Normalized', 'Position',[0.75 0.07 0.2 0.08 ], 'String','Remove Selected', 'CallBack', {@btnRemoveStandard_Callback});        
        if isfield(NMRdata, 'standards')
            set(NMRdata.loadStandardsData.lbStandardList, 'string', NMRdata.standards.label{1});
            plotSpectra('StdFigure', 1);
        end
        if ~isfield(NMRdata, 'standards') || isempty(NMRdata.standards.data),
            set(NMRdata.loadStandardsData.btnRemoveStandard,'Enable','off');
        end
        set(NMRdata.hLoadStandards,  'Visible', 'on');               
    end

%% UIcontrol callbacks: Options window
    function btnSaveGeneralOptions_Callback (~,~)
        NMRdata.fitPars.maxError=get(NMRdata.generalOptionsData.txtMaxError, 'string');
        str=['##MAXERROR=' NMRdata.fitPars.maxError  '\n'];
        NMRdata.fitPars.maxError=str2double(NMRdata.fitPars.maxError);
        NMRdata.fitPars.maxIters=get(NMRdata.generalOptionsData.txtmaxIters, 'string');
        str=[str '##MAXITERS=' NMRdata.fitPars.maxIters  '\n'];
        NMRdata.fitPars.maxIters=str2double(NMRdata.fitPars.maxIters);
        NMRdata.fitPars.normAcqPars=get(NMRdata.generalOptionsData.cbNormACQPars, 'value');
        str=[str '##NORMACQPARS=' num2str(NMRdata.fitPars.normAcqPars)  '\n'];
        NMRdata.paths.patterns=get(NMRdata.generalOptionsData.txtPathPatterns, 'string');
        str=[str '##PATTERNS_FOLDER=' strrep(NMRdata.paths.patterns, '\', '\\')  '\n'];
        NMRdata.paths.stds=get(NMRdata.generalOptionsData.txtPathStd, 'string');
        str=[str '##STDS_FOLDER=' strrep(NMRdata.paths.stds, '\', '\\')  '\n'];
        if isfield(NMRdata, 'standards'),
            NMRdata=rmfield(NMRdata, 'standards');
        end        
        LoadFileStandard(NMRdata.paths.stds, []);
        f = fopen('options.nmrcfg', 'w+');
        fprintf(f, str);
        fclose(f);        
        close(NMRdata.hEditGeneralOptions);
    end

%% UIcontrol callbacks: Main window
    function lbVisible_Callback (hObject,~)
        visibleIndex=get(hObject, 'value');
        if visibleIndex(1)==1,
            set(hObject, 'value', 1);
            NMRdata.visible=1:size(NMRdata.X.data, 1);
        else
            NMRdata.visible=visibleIndex-1;
        end
        plotSpectra('MainFigure', NMRdata.visible);
    end

    function lblMainXaxis_Callback (~,~)
        str=get(NMRdata.lblMainXaxis, 'string');
        if strcmp(str, 'chemical shift (ppm)')
            NMRdata.scaleXaxis{2}=NMRdata.scaleXaxis{2}*NMRdata.X.SF(1)/1e6;
            set(NMRdata.lblMainXaxis, 'string', 'frequency (Hz)');
        else
            NMRdata.scaleXaxis{2}=NMRdata.X.axisscale{2};
            set(NMRdata.lblMainXaxis, 'string', 'chemical shift (ppm)');
        end
        plotSpectra('MainFigure', NMRdata.visible);
    end

    function tblIncludedSamples_CellEditCallback (hObject,~)
        data=get(hObject, 'Data');
        NMRdata.X.include{1}=find(cell2mat(data(:,1)));
    end

    function btnSelAll_Callback (~,~)
        old=get(NMRdata.uiControls.tblIncludedSamples, 'data');
        if all(cell2mat(old(:,1))),
            included=num2cell(false(size(NMRdata.X.data,1),1));
            NMRdata.X.include{1}=false(size(NMRdata.X.data,1),1);
        else
            included=num2cell(true(size(NMRdata.X.data,1),1));
            NMRdata.X.include{1}=find(true(size(NMRdata.X.data,1),1));
        end
        cellArray=[included, NMRdata.X.label{1}];
        set(NMRdata.uiControls.tblIncludedSamples, 'data', cellArray);        
    end

    function popRefSample_Callback (hObject,~)
        NMRdata.fitPars.correctShape.sampleRef=get(hObject, 'value');
    end

    function popCorrectWidth_Callback (hObject, ~)
        aux=get(hObject, 'value');
        listStr=get(hObject, 'string');
        RefMode=listStr{aux};
        NMRdata.fitPars.correctShape.widthRef=RefMode;        
        if strcmp('None', RefMode) || ~isfield(NMRdata, 'X'),
            set(eventdata, 'value', 1);
            set(NMRdata.uiControls.popRefSample, 'value', 1);
            NMRdata.fitPars.correctShape.sampleRef=1;
            set(NMRdata.uiControls.lblRefSample, 'visible', 'off');
            set(NMRdata.uiControls.popRefSample, 'visible', 'off');
        else
            set(NMRdata.uiControls.lblRefSample, 'visible', 'on');
            set(NMRdata.uiControls.popRefSample, 'string', NMRdata.X.label{1});
            set(NMRdata.uiControls.popRefSample, 'visible', 'on');
            if strcmp('off', get(NMRdata.uiControls.popRefSample, 'visible'))
                set(NMRdata.uiControls.popRefSample, 'value', 1);
                NMRdata.fitPars.correctShape.sampleRef=1;
            end            
        end
    end

    function btnFit_Callback (~,~)
        bCorrectShape=1;
        nbSamples=size(NMRdata.X.data,1);
        nbIncSamples=length(NMRdata.X.include{1}==1);
        IncCount=0;
        if any(NMRdata.X.include{1})
            correctWidthFactor=ones(1,nbSamples);
            switch NMRdata.fitPars.correctShape.widthRef
                case 'None'
                    bCorrectShape=0;
                otherwise
                    bCorrectShape=0;
            end
            if bCorrectShape==1,            
                indexFWHH=strcmp(NMRdata.parResults, 'Width at half-height (Hz)');
                resultsFWHH=get(NMRdata.uiControls.tblresults(indexFWHH), 'data');
                signal2CorrectFWHH=get(NMRdata.uiControls.popCorrectWidth, 'value')-1; 
                correctWidthFactor=resultsFWHH(:, signal2CorrectFWHH)/resultsFWHH(NMRdata.fitPars.correctShape.sampleRef, signal2CorrectFWHH);
            end        
            nbPatterns=length(NMRdata.pattern);
            breakAll=0;
            wb = waitbar(0/nbIncSamples, ['Fitting region: ' '0/' num2str(nbPatterns) '. Please wait...'], 'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
            setappdata(wb,'canceling',0)
            for l=1:nbPatterns,
                currentPattern=NMRdata.pattern(l);
                currentPattern.mets=[];
                protons=[];
                range=getRang(NMRdata.X.axisscale{2}, NMRdata.pattern(l).maxppm, NMRdata.pattern(l).minppm);
                parseOK=1;
                if NMRdata.pattern(l).nbMets==0,
                    parseOK=0;
                elseif currentPattern.qMode==2,
                    if NMRdata.pattern(l).mets.met1.type==1,
                        tag='met1';
                    else
                        tag='spec1';
                    end
                    currentPattern.nbMets=1;
                    currentPattern.mets.(tag).name=NMRdata.pattern(l).name;
                    protons=NMRdata.pattern(l).mets.met1.protons;
                else
                    for j=1:NMRdata.pattern(l).nbMets,
                        protons(j)=NMRdata.pattern(l).mets.(['met' num2str(j)]).protons;
                        if NMRdata.pattern(l).mets.(['met' num2str(j)]).type==1,    
                            currentPattern.mets.(['met' num2str(j)])=NMRdata.pattern(l).mets.(['met' num2str(j)]);
                        else
                            currentPattern.mets.(['spec' num2str(j)])=NMRdata.pattern(l).mets.(['met' num2str(j)]);
                            stdIndex=strcmp(currentPattern.mets.(['spec' num2str(j)]).standard, NMRdata.standards.label{1});
                            if ~isempty(stdIndex) && any(stdIndex),
                                rangeSpec=getRang(NMRdata.X.axisscale{2}, currentPattern.mets.(['spec' num2str(j)]).hlim, currentPattern.mets.(['spec' num2str(j)]).llim);
                                currentPattern.mets.(['spec' num2str(j)]).xaxis=NMRdata.X.axisscale{2}(rangeSpec);
                                currentPattern.mets.(['spec' num2str(j)]).data=NMRdata.standards.data(stdIndex, rangeSpec);
                                %Subtract possible offset
                                currentPattern.mets.(['spec' num2str(j)]).data=currentPattern.mets.(['spec' num2str(j)]).data-min(currentPattern.mets.(['spec' num2str(j)]).data);
                                %normalise intensities
                                currentPattern.mets.(['spec' num2str(j)]).scalefactor=max(currentPattern.mets.(['spec' num2str(j)]).data);
                                currentPattern.mets.(['spec' num2str(j)]).data=currentPattern.mets.(['spec' num2str(j)]).data/currentPattern.mets.(['spec' num2str(j)]).scalefactor;
                            else
                                warndlg(['The standard "'  currentPattern.mets.(['spec' num2str(j)]).standard{1} '" for signal "' currentPattern.mets.(['spec' num2str(j)]).name '" was not found. You must load it to current session before fitting this pattern']);
                                parseOK=0;
                            end
                        end
                    end
                end
                if parseOK==1,
                    results=zeros(nbSamples, length(fieldnames(currentPattern.mets)), 6);
                    waitbar(0/nbIncSamples, wb, ['Fitting region: ' num2str(l) '/' num2str(nbPatterns) '. Please wait...']);
                    IncCount=0;
                    for j=1:nbSamples,
                        if getappdata(wb,'canceling'),
                            delete(wb);
                            breakAll=1;
                            break;
                        end
                        waitbar(IncCount/nbIncSamples, wb, ['Fitting region: ' num2str(l) '/' num2str(nbPatterns) '. Please wait...']);
                        if ismember(j, NMRdata.X.include{1}),
                            for k=1:currentPattern.nbMets,
                                %RbR: Unificar
                                if NMRdata.pattern(l).mets.(['met' num2str(k)]).type==1,
                                    currentPattern.mets.(['met' num2str(k)]).fwhh=NMRdata.pattern(l).mets.(['met' num2str(k)]).fwhh*correctWidthFactor(j);
                                elseif NMRdata.pattern(l).mets.(['met' num2str(k)]).type==2,
                                    currentPattern.mets.(['spec' num2str(k)]).fwhh=NMRdata.pattern(l).mets.(['met' num2str(k)]).fwhh*correctWidthFactor(j);
                                end
                            end
                            plotMode=0;
                            plotsPath=get(NMRdata.uiControls.txtFigPath, 'String');
                            if exist(plotsPath,'dir'),
                                plotMode=1;
                            end
                            x.data=NMRdata.X.data(j, :);
                            x.axisscale{2}=NMRdata.X.axisscale{2};
                            adParams.xaxisscale=NMRdata.X.axisscale{2};
                            xNoise = changeSpectralRange( x, [x.axisscale{2}(1) (x.axisscale{2}(1)-0.5)], 'i', 2 );
                            noiseLevel = getNoiseLevel(xNoise.data(1,:));
                            x=changeSpectralRange( x, [NMRdata.pattern(l).maxppm NMRdata.pattern(l).minppm], 'i', 2 );
                            tic
                            if currentPattern.qMode==1,
                                try
                                    [ output_args ] = fitNMRsignals( x.data, x.axisscale{2}, noiseLevel, currentPattern.mets, NMRdata.X.SF(1), currentPattern.baseline, currentPattern.constraints', NMRdata.fitPars.maxError, NMRdata.fitPars.maxIters, 0);
                                catch
                                    delete(wb);
                                end
                            else
                                output_args.pureSpectra=sum(sum(real(x.data)));
                                output_args.h=figure;
                                set(output_args.h, 'visible', 'off');
                                output_args.s(1:5,:)=0;
                                plot(x.axisscale{2}, real(x.data));
                                protons=1;
                            end
                            toc
                            if plotMode==1,
                                print(output_args.h, [plotsPath '\' NMRdata.X.label{1}{j,:} '_' currentPattern.name], '-dpng', '-r300');
                                close(output_args.h);
                            end
                            results(j, :, 2:end)=output_args.s(1:5,:)';  
                            results(j, :, 1)=sum(sum(output_args.pureSpectra,3),2)./protons';
                            IncCount=IncCount+1;
                        end
                    end
                    if exist('wb')
                        if breakAll,
                            break;
                        elseif l==nbPatterns,
                            delete(wb);
                        end
                    end
                    submets=fieldnames(currentPattern.mets);
                    n1=length(submets);
                    n2=length(NMRdata.uiControls.tblresults);
                    colHeaders=get(NMRdata.uiControls.tblresults(1), 'columnname');
                    resultMatrix=get(NMRdata.uiControls.tblresults(1), 'data');
                    [r c]=size(resultMatrix);
                    resultMatrix=zeros(nbSamples, c, length(NMRdata.uiControls.tblresults));
                    indexMatrix=zeros(nbSamples, c);
                    for k=1:n1,
                        matchIndex=find(strcmp(colHeaders, {currentPattern.mets.(submets{k}).name}), 1);
                        if isempty(matchIndex)
                            colHeaders=[colHeaders; {currentPattern.mets.(submets{k}).name}];
                            matchIndex=length(colHeaders);
                            resultMatrix(1:nbSamples, matchIndex, 1:n2)=nan;
                        end
                        resultMatrix(NMRdata.X.include{1}, matchIndex, :)=results(NMRdata.X.include{1}, k, :);
                        indexMatrix(NMRdata.X.include{1}, matchIndex)=1;
                    end
                    indexMatrix=~indexMatrix;
                    for k=1:length(NMRdata.uiControls.tblresults),                
                        newData=squeeze(resultMatrix(:,:,k));
                        oldData=get(NMRdata.uiControls.tblresults(k), 'data');
                        nbElements=numel(oldData);
                        linearIndexMatrix=find(indexMatrix);
                        linearIndexMatrix(linearIndexMatrix>nbElements)=[];
                        newData(linearIndexMatrix)=oldData(linearIndexMatrix);
                        set(NMRdata.uiControls.tblresults(k), 'data', newData);
                        set(NMRdata.uiControls.tblresults(k), 'columnname', colHeaders);
                    end
                    widthCorrSelected=get(NMRdata.uiControls.popCorrectWidth, 'value');           
                    s=get(NMRdata.uiControls.popCorrectWidth, 'String');
                    widthCorrSelected=find(strcmp(s(widthCorrSelected), [{'None'}; colHeaders]));
                    set(NMRdata.uiControls.popCorrectWidth, 'String', [{'None'}; colHeaders]);
                    set(NMRdata.uiControls.popCorrectWidth, 'value', widthCorrSelected);   
                else
                    delete(wb);
                end
            end
            if parseOK==1 && ~isempty(get(NMRdata.uiControls.tblresults(k), 'data')),
                NMRdata.parResults=[{'area'} {'intensity'} {'center (ppm)'} {'Width at half-height (Hz)'} {'Gaussian contribution'} {'J-coupling (Hz)'} ];
                set(NMRdata.uiControls.lblParResult, 'string', NMRdata.parResults(end));
                guidata(hMainFigure, NMRdata);
                btnParResult_Callback(0,0,1);
            end
        end
    end

%Pattern settings
    function btnLoadPattern_Callback (~,~)
        [filename, pathname] = uigetfile([NMRdata.paths.patterns '\*.nmrsgnl']);
        if filename<1
            return;
        end
        f = fopen([pathname filename]); 
        patternFile = fread(f);
        fclose(f);
        patternFile = char(patternFile');
        patternFile=regexprep(patternFile,'\r\n|\n|\r',''); 
        if ~strcmp(patternFile,'')          
            posCR=[strfind(patternFile, sprintf('##')) length(patternFile)+1];
            line=1;
            iPattern=1;
            if isfield(NMRdata, 'pattern')
                iPattern=length(NMRdata.pattern)+1;
            end
            mode=3;
            standardPointer=1;
            while (line<length(posCR))
                currentLine=patternFile(posCR(line):posCR(line+1)-1);
                if strcmp(currentLine(1:2), '##')
                    aux=strfind(currentLine, '=');
                    switch currentLine(3:aux(1)-1)
                        case 'NAME'
                            NMRdata.pattern(iPattern).name=(currentLine(aux+1:end));
                        case 'MINPPM'
                            NMRdata.pattern(iPattern).minppm=str2double(currentLine(aux+1:end));
                        case 'MAXPPM'
                            NMRdata.pattern(iPattern).maxppm=str2double(currentLine(aux+1:end));
                        case 'BASEMODE'
                            NMRdata.pattern(iPattern).baseline=str2double(currentLine(aux+1:end));                        
                        case 'QMODE'
                            switch currentLine(aux+1:end)
                                case 'LS'
                                    mode=1;
                                case 'INT'
                                    mode=2;
                            end
                            NMRdata.pattern(iPattern).qMode=mode;
                        case 'nbMETS'
                            NMRdata.pattern(iPattern).nbMets=str2double(currentLine(aux+1:end));
                        case 'METS'
                            labels=textscan(currentLine(aux+1:end), '%s', 'Delimiter', ',');
                            labels=labels{1};
                            for i=1:NMRdata.pattern(iPattern).nbMets,
                                NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).name=labels{i};
                            end
                        case 'STANDARDS'
                            standardNames=textscan(currentLine(aux+1:end), '%s', 'Delimiter', ',');      
                            standardNames=standardNames{1};
                        case 'PARAMS'
                            aux=textscan(currentLine(aux+1:end), '%f', 'Delimiter', ',');
                            aux=aux{1};
                            aux=reshape(aux', 7, NMRdata.pattern(iPattern).nbMets)';
                            for i=1:NMRdata.pattern(iPattern).nbMets,
                                NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).type=aux(i,1);
                                NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).pos=aux(i,2);
                                NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).fwhh=aux(i,3);
                                if NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).type==1,                                
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).gauss=aux(i,4);
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).jcop=aux(i,5);                                    
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).mult=aux(i,6);
                                else
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).stdPos=aux(i,4);
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).hlim=aux(i,5);
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).llim=aux(i,6);
                                    NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).standard=standardNames(standardPointer);
                                    standardPointer=standardPointer+1;
                                end
                                NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).protons=aux(i,7);
                            end
                        case 'CONSTRAINTS'
                            aux=textscan(currentLine(aux+1:end), '%f', 'Delimiter', ',');
                            aux=aux{1};
                            aux=reshape(aux', 5, NMRdata.pattern(iPattern).nbMets)';
                            NMRdata.pattern(iPattern).constraints=aux;
                    end
                end
                line=line+1;
            end
            regions=extractfield(NMRdata.pattern, 'name');
            set(NMRdata.uiControls.lbPatterns, 'string', regions);
            set(NMRdata.uiControls.lbPatterns, 'value', iPattern);
            lbPatterns_Callback(0,0,iPattern);
            set(NMRdata.uiControls.btnDelPattern, 'enable', 'on');
            set(NMRdata.uiControls.btnSavePattern, 'enable', 'on');
            h=get(NMRdata.pnlRegionSettings, 'children');
            EnableDisableUIcontrols(h, 'on');
        end
    end

    function btnSavePattern_Callback (~,~)
        fitmode=[{'LS'},{'INT'}];
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        baseline=NMRdata.pattern(iPattern).baseline;
        nbMets=NMRdata.pattern(iPattern).nbMets;
        strStandards='';
        aux=zeros(nbMets, 7);
        name=NMRdata.pattern(iPattern).name;
        str=['##NAME=' name  '\n'];
        str=[str '##MINPPM=' num2str(NMRdata.pattern(iPattern).minppm) '\n'];
        str=[str '##MAXPPM=' num2str(NMRdata.pattern(iPattern).maxppm) '\n'];
        str=[str '##QMODE=' fitmode{NMRdata.pattern(iPattern).qMode} '\n'];
        if strcmp(fitmode{NMRdata.pattern(iPattern).qMode}, 'INT'),
            baseline=0;
            nbMets=0;
            constraints=0;
            aux=0;
            
        end
        str=[str '##BASEMODE=' num2str(baseline) '\n'];
        str=[str '##nbMETS=' num2str(nbMets) '\n'];
        names=[];
        for i=1:nbMets,
            names=[names ',' NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).name];
            aux(i,1)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).type;
            aux(i,2)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).pos;
            aux(i,3)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).fwhh;
            if NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).type==1,
                aux(i,4)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).gauss;
                aux(i,5)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).jcop;
                aux(i,6)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).mult;
            else
                strStandards=[strStandards ',' char(NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).standard)];
                aux(i,4)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).stdPos;
                aux(i,5)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).hlim;
                aux(i,6)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).llim;
            end
            aux(i,7)=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).protons;
        end
        if nbMets>0,
            constraints=NMRdata.pattern(iPattern).constraints;            
            str=[str '##METS=' names(2:end) '\n'];
            if ~isempty(strStandards),
                str=[str '##STANDARDS=' strStandards(2:end) '\n'];
            end
            str=[str '##PARAMS=\n'];
            aux=mat2str(aux);
            aux=strrep(aux, ';', ',\n');
            aux=strrep(aux, ' ', ',');
            str=[str aux(2:end-1) '\n'];
            str=[str '##CONSTRAINTS=\n'];
            aux=mat2str(constraints);
            aux=strrep(aux, ';', ',\n');
            aux=strrep(aux, ' ', ',');
            str=[str aux(2:end-1) '\n'];   
        end
        [fileName, filepath] = uiputfile('*.nmrsgnl', 'Save As', [NMRdata.paths.patterns '\' NMRdata.pattern(iPattern).name]);
        f = fopen([filepath fileName], 'w+');
        fprintf(f, str);
        fclose(f);
    end

    function btnDelPattern_Callback(~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        if ~isempty(iPattern)
            list=get(NMRdata.uiControls.lbPatterns, 'string');
            list(iPattern)=[];
            set(NMRdata.uiControls.lbPatterns, 'string', list);
            set(NMRdata.uiControls.lbPatterns, 'value', 1);            
            if length(NMRdata.pattern)==1,
                NMRdata=rmfield(NMRdata, 'pattern');   
                iPattern=0;
            else
                NMRdata.pattern(iPattern)=[];
                iPattern=1;
            end
        end
        lbPatterns_Callback(0,0,iPattern);  
    end

    function btnCreatePattern_Callback(~,~)
        name=get(NMRdata.uiControls.txtCreateRegion, 'string');
        set(NMRdata.uiControls.txtCreateRegion, 'string', '');
        if ~isempty(name),
            iPattern=1;
            if isfield(NMRdata, 'pattern')
                iPattern=iPattern+length(NMRdata.pattern);  
            end
            NMRdata.pattern(iPattern).name=name;
            NMRdata.pattern(iPattern).minppm=0;
            NMRdata.pattern(iPattern).maxppm=0;
            NMRdata.pattern(iPattern).baseline=0;
            NMRdata.pattern(iPattern).qMode=1;
            NMRdata.pattern(iPattern).nbMets=0;
            list=get(NMRdata.uiControls.lbPatterns, 'string');
            list=[list; {name}];
            set(NMRdata.uiControls.lbPatterns, 'string', list); 
            set(NMRdata.uiControls.lbPatterns, 'value', iPattern);
            guidata(hMainFigure, NMRdata); 
            lbPatterns_Callback(0, 0, iPattern);
            set(NMRdata.uiControls.btnDelPattern, 'enable', 'on');
            set(NMRdata.uiControls.btnSavePattern, 'enable', 'on');            
            h=get(NMRdata.pnlRegionSettings, 'children');
            EnableDisableUIcontrols(h, 'on');                        
        end 
    end

    function lbPatterns_Callback (hObject, ~, userData)
        if isfield(NMRdata,'hEditMet')
            close(NMRdata.hEditMet)
        end        
        if hObject~=0,
            iPattern=get(hObject, 'value');
        else
            iPattern=userData;
        end
        if iPattern>0
            set(NMRdata.uiControls.txtMaxPPM, 'string', num2str(NMRdata.pattern(iPattern).maxppm));
            set(NMRdata.uiControls.txtMinPPM, 'string', num2str(NMRdata.pattern(iPattern).minppm));
            if bitand(1, NMRdata.pattern(iPattern).baseline),
                set(NMRdata.uiControls.cbFitBaselinePol, 'value', 1);
            else
                set(NMRdata.uiControls.cbFitBaselinePol, 'value', 0);
            end
            if bitand(2, NMRdata.pattern(iPattern).baseline),
                set(NMRdata.uiControls.cbFitBaselineGauss, 'value', 1);
            else
                set(NMRdata.uiControls.cbFitBaselineGauss, 'value', 0);
            end
            if bitand(4, NMRdata.pattern(iPattern).baseline),
                set(NMRdata.uiControls.cbFitBaselineCos, 'value', 1);
            else
                set(NMRdata.uiControls.cbFitBaselineCos, 'value', 0);
            end            
            set(NMRdata.uiControls.popFitMode, 'value', NMRdata.pattern(iPattern).qMode);
            nameMets=[];
            for i=1:NMRdata.pattern(iPattern).nbMets,
                nameMets{i}=NMRdata.pattern(iPattern).mets.(['met' num2str(i)]).name;
            end
            set(NMRdata.uiControls.lbMetList, 'string', nameMets);
            set(NMRdata.uiControls.lbMetList, 'value', 1);
        else
            set(NMRdata.uiControls.txtMaxPPM, 'string', '');
            set(NMRdata.uiControls.txtMinPPM, 'string', '');
            set(NMRdata.uiControls.cbFitBaselinePol, 'value', 0);
            set(NMRdata.uiControls.cbFitBaselineGauss, 'value', 0);
            set(NMRdata.uiControls.cbFitBaselineCos, 'value', 0);
            set(NMRdata.uiControls.lbMetList, 'string', []);
            set(NMRdata.uiControls.btnDelPattern, 'enable', 'off');
            set(NMRdata.uiControls.btnSavePattern, 'enable', 'off');            
            h=get(NMRdata.pnlRegionSettings, 'children');
            EnableDisableUIcontrols(h, 'off');
        end
    end

    function popFitMode_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        NMRdata.pattern(iPattern).qMode=get(NMRdata.uiControls.popFitMode, 'value');
    end

    function cbFitBaseline_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        bitBase=0;
        bitBase=bitor(bitBase, get(NMRdata.uiControls.cbFitBaselinePol, 'value')*1);
        bitBase=bitor(bitBase, get(NMRdata.uiControls.cbFitBaselineGauss, 'value')*2);
        bitBase=bitor(bitBase, get(NMRdata.uiControls.cbFitBaselineCos, 'value')*4);
        NMRdata.pattern(iPattern).baseline=bitBase;
    end

    function lblMaxPPM_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        NMRdata.pattern(iPattern).maxppm=str2double(get(NMRdata.uiControls.txtMaxPPM, 'String'));
    end

    function lblMinPPM_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        NMRdata.pattern(iPattern).minppm=str2double(get(NMRdata.uiControls.txtMinPPM, 'String'));
    end

%Signal settings
    function sldConstraint_Callback (hObject,~)
        value=get(hObject,'value');
        if rem(value,1)~=0,
            value=round(value);
            set(hObject, 'value', value);
        end
        switch hObject
            case NMRdata.editMetData.sldConstCenter
                factor=0.001;
                handle=NMRdata.editMetData.lblConstCenter;
            case NMRdata.editMetData.sldConstFWHH
                factor=0.05;
                handle=NMRdata.editMetData.lblConstFWHH;
            case NMRdata.editMetData.sldConstGaussian
                factor=0.1;
                handle=NMRdata.editMetData.lblConstGaussian;                
            case NMRdata.editMetData.sldConstJcoupling
                factor=0.1;
                handle=NMRdata.editMetData.lblConstJcoupling;
        end
        value=value*factor;
        value=num2str(value);
        set(handle, 'string', [char(177) value]);
    end

    function btnMetSave_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        value=str2double(get(NMRdata.editMetData.txtCenter, 'String'));
        NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).pos=value;
        value=get(NMRdata.editMetData.lblConstCenter, 'String');
        value=str2double(value(2:end))/0.001;
        NMRdata.pattern(iPattern).constraints(NMRdata.editMetData.metIndex, 2)=value;
        value=str2double(get(NMRdata.editMetData.txtFWHH, 'String'));
        NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).fwhh=value;
        value=get(NMRdata.editMetData.lblConstFWHH, 'String');
        value=str2double(value(2:end))/0.05;        
        NMRdata.pattern(iPattern).constraints(NMRdata.editMetData.metIndex, 3)=value;
        value=get(NMRdata.editMetData.lblConstGaussian, 'String');
        value=str2double(value(2:end))/0.1;            
        NMRdata.pattern(iPattern).constraints(NMRdata.editMetData.metIndex, 4)=value;      
        if NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).type==1,
            value=str2double(get(NMRdata.editMetData.txtGaussian, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).gauss=value;         
            value=str2double(get(NMRdata.editMetData.txtJcoupling, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).jcop=value;   
            value=get(NMRdata.editMetData.lblConstJcoupling, 'String');
            value=str2double(value(2:end))/0.1;           
            NMRdata.pattern(iPattern).constraints(NMRdata.editMetData.metIndex, 5)=value; 
            if get(NMRdata.editMetData.cbApplyAllMets, 'value'),
                NMRdata.pattern(iPattern).constraints=repmat(NMRdata.pattern(iPattern).constraints(NMRdata.editMetData.metIndex,:), NMRdata.pattern(iPattern).nbMets, 1);
            end
            value=str2double(get(NMRdata.editMetData.txtMultiplicity, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).mult=value;            
        else     
            selStd=get(NMRdata.editMetData.popNameStd, 'value');
            selStdName=get(NMRdata.editMetData.popNameStd, 'string');
            selStdName=selStdName(selStd);
            if ~strcmp(selStdName, ' '),
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).standard=selStdName;
            end
            value=str2double(get(NMRdata.editMetData.txtCenterStd, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).stdPos=value;            
            value=str2double(get(NMRdata.editMetData.txtMaxLimitStd, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).hlim=value;
            value=str2double(get(NMRdata.editMetData.txtMinLimitStd, 'String'));
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).llim=value;            
        end
        value=str2double(get(NMRdata.editMetData.txtNbProtons, 'String'));
        NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.editMetData.metIndex)]).protons=value;            
        close(NMRdata.hEditMet);
    end

    function btnDelMet_Callback (~,~)
        iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
        index2del=get(NMRdata.uiControls.lbMetList, 'value');
        list=get(NMRdata.uiControls.lbMetList, 'string');
        for j=index2del:NMRdata.pattern(iPattern).nbMets-1,
            NMRdata.pattern(iPattern).mets.(['met' num2str(j)])=NMRdata.pattern(iPattern).mets.(['met' num2str(j+1)]);
        end
        NMRdata.pattern(iPattern).mets=rmfield(NMRdata.pattern(iPattern).mets, ['met' num2str(NMRdata.pattern(iPattern).nbMets)]);
        NMRdata.pattern(iPattern).constraints(index2del,:)=[];
        NMRdata.pattern(iPattern).nbMets=NMRdata.pattern(iPattern).nbMets-1;
        list(index2del)=[];
        set(NMRdata.uiControls.lbMetList, 'string', list);
        set(NMRdata.uiControls.lbMetList, 'value', 1);
    end

    function btnNewMet_Callback (~,~)
        name=get(NMRdata.uiControls.txtNewMet, 'string');
        list=get(NMRdata.uiControls.lbMetList, 'string');
        norepeated=1;
        if ~isempty(list),
            norepeated=isempty(find(cellfun(@(x) (strcmp(name, x)), list ), 1));
        end
        if ~isempty(name) && norepeated,
            type=questdlg('Choose signal type','New Signal', 'Voigt model','Standard template', '');
            switch type
                case 'Voigt model',
                    type=1;
                case 'Standard template',
                    if isfield(NMRdata, 'standards');
                        stdSel = listdlg('ListString', NMRdata.standards.label{1}, 'SelectionMode','single');
                        if isempty(stdSel)
                            return;
                        end
                        type=2;
                    else
                        warndlg('There are no standards currently loaded');
                        return;
                    end
                otherwise
                    type=-1;
            end
            if type<1
                return;
            end
            iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).type=type;
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).name=name;
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).pos=0;  
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).fwhh=1;
            if type==1,
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).jcop=0;
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).gauss=0;            
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).mult=1;
            else
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).standard=NMRdata.standards.label{1}(stdSel);
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).stdPos=0;
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).hlim=0;
                NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).llim=0;
            end
            NMRdata.pattern(iPattern).mets.(['met' num2str(NMRdata.pattern(iPattern).nbMets+1)]).protons=1;
            NMRdata.pattern(iPattern).constraints(NMRdata.pattern(iPattern).nbMets+1, 1:5)=0;
            list=[list; {name}];
            set(NMRdata.uiControls.lbMetList, 'string', list);
            NMRdata.pattern(iPattern).nbMets=NMRdata.pattern(iPattern).nbMets+1;
            set(NMRdata.uiControls.lbMetList, 'value', NMRdata.pattern(iPattern).nbMets);           
        end
        set(NMRdata.uiControls.txtNewMet, 'string', '');        
    end

    function lbMetList_Callback (hObject, ~)
        metIndex=get(hObject, 'value');
        if metIndex>0
            if isfield(NMRdata,'hEditMet')
                close(NMRdata.hEditMet)
            end
            NMRdata.editMetData.metIndex=metIndex;
            iPattern=get(NMRdata.uiControls.lbPatterns, 'value');
            NMRdata.hEditMet=figure('Units','pixels', 'NumberTitle', 'off', 'Name', 'Modify signal properties', 'MenuBar','none', 'Toolbar','none', 'DeleteFcn',@hEditMet_DeleteFcn, 'Visible','off');
            pnlOptionsPanel=uipanel(NMRdata.hEditMet, 'FontWeight','bold', 'TitlePosition','centertop', 'Units','Normalized', 'Position',[0.05 0.05 0.9 0.9]);
            lblInitValues=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'FontWeight','bold', 'HorizontalAlignment','center', 'Position',[0.225 0.8 0.15 0.05], 'String', 'Init values' );                 
            lblConstraints=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'FontWeight','bold', 'HorizontalAlignment','center', 'Position',[0.675 0.8 0.15 0.05], 'String', 'Constraints' );                 
            lblCenter=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.7 0.28 0.05], 'String', 'Center (ppm): ' );         
            NMRdata.editMetData.txtCenter=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.7 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).pos));
            currentValue=NMRdata.pattern(iPattern).constraints( metIndex,2);
            NMRdata.editMetData.sldConstCenter = uicontrol(pnlOptionsPanel, 'Style','slider', 'Units','Normalized', 'Min',0, 'Max',100, 'Value',currentValue, 'SliderStep',[0.01 0.1], 'Position', [0.6 0.7 0.25 0.05],  'CallBack', {@sldConstraint_Callback});
            NMRdata.editMetData.lblConstCenter = uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','center', 'Position', [0.85 0.7 0.1 0.05], 'String',[char(177) num2str(currentValue*0.001)]);
            NMRdata.editMetData.txtFWHH=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.6 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).fwhh));
            currentValue=NMRdata.pattern(iPattern).constraints( metIndex,3);        
            NMRdata.editMetData.sldConstFWHH= uicontrol(pnlOptionsPanel, 'Style', 'slider', 'Units','Normalized', 'Min',0,'Max',50,'Value',currentValue, 'SliderStep',[0.02 0.2], 'Position', [0.6 0.6 0.25 0.05], 'CallBack', {@sldConstraint_Callback});        
            NMRdata.editMetData.lblConstFWHH= uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','center', 'Position', [0.85 0.6 0.1 0.05], 'String',[char(177) num2str(currentValue*0.05)]);        
            lblGaussian=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.5 0.3 0.05], 'String', 'Gaussian contribution (0-1): ' );                 
            currentValue=NMRdata.pattern(iPattern).constraints( metIndex,4);         
            NMRdata.editMetData.sldConstGaussian = uicontrol(pnlOptionsPanel, 'Style', 'slider', 'Units','Normalized', 'Min',0,'Max',10,'Value',currentValue, 'SliderStep',[0.1 0.2], 'Position', [0.6 0.5 0.25 0.05], 'CallBack', {@sldConstraint_Callback});        
            NMRdata.editMetData.lblConstGaussian = uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','center', 'Position', [0.85 0.5 0.1 0.05], 'String',[char(177) num2str(currentValue*0.1)]);                        
            if NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).type==1,
                lblFWHH=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.6 0.28 0.05], 'String', 'Width at half-height (Hz): ' );         
                NMRdata.editMetData.txtGaussian=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.5 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).gauss));
                lblJcoupling=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.4 0.28 0.05], 'String', 'J-coupling (Hz): ' );         
                NMRdata.editMetData.txtJcoupling=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.4 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).jcop));
                currentValue=NMRdata.pattern(iPattern).constraints( metIndex,5);                
                NMRdata.editMetData.sldConstJcoupling = uicontrol(pnlOptionsPanel, 'Style', 'slider', 'Units','Normalized', 'Min',0,'Max',10,'Value',currentValue, 'SliderStep',[0.1 0.2], 'Position', [0.6 0.4 0.25 0.05], 'CallBack', {@sldConstraint_Callback});        
                NMRdata.editMetData.lblConstJcoupling = uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','center', 'Position', [0.85 0.4 0.1 0.05], 'String',[char(177) num2str(currentValue*0.1)]);        
                lblMultiplicity=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.3 0.28 0.05], 'String', 'Multiplicity: ' );         
                NMRdata.editMetData.txtMultiplicity=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.3 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).mult));                               
                lblNbProtons=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.2 0.28 0.05], 'String', 'Number of protons: ' );         
                NMRdata.editMetData.txtNbProtons=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.2 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).protons));        
                lblApplyAllMets=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.6 0.25 0.3 0.1], 'String','Copy these constraints to all signals: ' ); 
                NMRdata.editMetData.cbApplyAllMets=uicontrol('Parent', pnlOptionsPanel, 'Style','checkbox', 'Units','Normalized', 'Position',[0.9 0.27 0.04 0.04]);                
            else
                lblFWHH=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.60 0.28 0.05], 'String', 'FWHH increase (Hz): ' );                         
                lblNameStd=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.9 0.28 0.05], 'String', 'Standard: ' ); 
                matchStd=0;
                listStds=[];
                if isfield(NMRdata, 'standards'),
                    matchStd=strcmp(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).standard{:}, NMRdata.standards.label{1});
                    listStds=NMRdata.standards.label{1};
                end
                if any(matchStd),
                    selStd=find(matchStd);
                else
                    listStds=[{' '}, listStds];
                    selStd=1;
                end
                NMRdata.editMetData.popNameStd=uicontrol('Parent',pnlOptionsPanel, 'Style','popup', 'Units','Normalized', 'backgroundColor', 'w', 'Position',[0.6 0.9 0.25 0.05], 'String', listStds, 'value', selStd);              
                lblCenterStd=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.4 0.28 0.05], 'String', 'standard center (ppm): ' );         
                NMRdata.editMetData.txtCenterStd=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.4 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).stdPos));                                                        
                lblMaxLimitStd=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.3 0.28 0.05], 'String', 'standard upper limit (ppm): ' );         
                NMRdata.editMetData.txtMaxLimitStd=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.3 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).hlim));                    
                lblMinLimitStd=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.2 0.28 0.05], 'String', 'standard lower limit (ppm): ' );         
                NMRdata.editMetData.txtMinLimitStd=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.2 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).llim));                                
                lblNbProtons=uicontrol(pnlOptionsPanel, 'Style','text', 'Units','Normalized', 'HorizontalAlignment','left', 'Position',[0.05 0.1 0.28 0.05], 'String', 'Number of protons: ' );         
                NMRdata.editMetData.txtNbProtons=uicontrol(pnlOptionsPanel, 'Style','edit', 'Units','Normalized', 'backgroundColor','w', 'Position',[0.35 0.1 0.1 0.05], 'String',num2str(NMRdata.pattern(iPattern).mets.(['met' num2str(metIndex)]).protons));                    
            end        
            btnMetSave=uicontrol(pnlOptionsPanel, 'Style','PushButton', 'Units','Normalized', 'FontWeight','bold', 'Position',[0.75 0.05 0.2 0.1 ], 'String', 'Save', 'CallBack', {@btnMetSave_Callback});                
            set(NMRdata.hEditMet,  'Visible', 'on');
        end
    end

%Result
    function tblresults_CellSelectionCallback (~,eventdata)
        NMRdata.tableColSel=eventdata.Indices(:,2);
    end

    function btnResultDelCol_Callback (~,~)
        if ~isempty(NMRdata.tableColSel)
            auxCol=get(NMRdata.uiControls.tblresults(1), 'columnname');
            auxCol(NMRdata.tableColSel)=[];
            for i=1:6,
                auxData=get(NMRdata.uiControls.tblresults(i), 'data');
                auxData(:, NMRdata.tableColSel)=[];
                set(NMRdata.uiControls.tblresults(i), 'data', auxData);
                set(NMRdata.uiControls.tblresults(i), 'columnname', auxCol);
                widthCorrSelected=get(NMRdata.uiControls.popCorrectWidth, 'value'); 
                s=get(NMRdata.uiControls.popCorrectWidth, 'String');
                widthCorrSelected=find(strcmp(s(widthCorrSelected), [{'None'}; auxCol]));
                if isempty(widthCorrSelected)
                    widthCorrSelected=1;
                end
                set(NMRdata.uiControls.popCorrectWidth, 'String', [{'None'}; auxCol]);
                set(NMRdata.uiControls.popCorrectWidth, 'value', widthCorrSelected);                
            end
        end
    end

    function btnResultExport_Callback (~,~)
        selPar=get(NMRdata.uiControls.lblParResult, 'string');
        name=datestr(now);
        name=strrep(name, ':', '_');
        name=[char(selPar) '_' name];        
        selPar=strcmp(NMRdata.parResults, selPar)';
        values=get(NMRdata.uiControls.tblresults(selPar), 'data');
        rows=get(NMRdata.uiControls.tblresults(selPar), 'rowname');
        columns=get(NMRdata.uiControls.tblresults(selPar), 'columnname');
        columns=[' ', columns'];
        [fileName, filepath] = uiputfile('*.csv', 'Save file', name);
        f = fopen([filepath fileName], 'w') ;
        fprintf(f, '%s\n', 'sep=;');
        fprintf(f, '%s;', columns{1:end});
        fprintf(f, '\n');
        for j=1:length(rows),
            fprintf(f, '%s;', rows{j});
            fclose(f);
            dlmwrite([filepath fileName], values(j,:), '-append', 'delimiter', ';');
            f = fopen([filepath fileName], 'r+') ;
            fseek(f, 0, 'eof');
        end
        fclose(f);
    end

    function btnParResult_Callback (~,~, increment)
        selPar=get(NMRdata.uiControls.lblParResult, 'string');
        selPar=strcmp(NMRdata.parResults, selPar)';
        if any(selPar)
            for i=1:length(NMRdata.uiControls.tblresults),
                set(NMRdata.uiControls.tblresults(i), 'visible', 'off');
            end
            selPar=circshift(selPar, [increment 0]);
            set(NMRdata.uiControls.tblresults(selPar), 'visible', 'on');
            set(NMRdata.uiControls.lblParResult, 'string', NMRdata.parResults(selPar));
        end
    end

%% UIcontrol callbacks: Standards' window
    function btnRemoveStandard_Callback (~,~)
        iStd2Rmv=get(NMRdata.loadStandardsData.lbStandardList, 'value');
        NMRdata.standards.label{1}(iStd2Rmv)=[];
        NMRdata.standards.data(iStd2Rmv,:)=[];
        set(NMRdata.loadStandardsData.lbStandardList, 'value', []);
        if isempty(NMRdata.standards.data),
            set(NMRdata.loadStandardsData.btnRemoveStandard,'Enable','off');
        end
        set(NMRdata.loadStandardsData.lbStandardList, 'string', NMRdata.standards.label{1});
        set(NMRdata.loadStandardsData.lbStandardList, 'value', 1);            
        plotSpectra('StdFigure', 1);
    end

    function btnLoadFileStandard_Callback (~,~)
        [filename, pathname] = uigetfile('*.*');
        if filename<1
            return;
        end 
        if ~strcmp(filename(end-6:end), '.nmrstd'),
            errordlg('Invalid File type');
            return;
        end
        LoadFileStandard(pathname, filename)           
        set(NMRdata.loadStandardsData.lbStandardList, 'string', NMRdata.standards.label{1});
        set(NMRdata.loadStandardsData.lbStandardList, 'value', 1);
        plotSpectra('StdFigure', 1);
        set(NMRdata.loadStandardsData.btnRemoveStandard,'Enable','on');
    end

    function lbStandardList_Callback (~,~)
        iStd2Show=get(NMRdata.loadStandardsData.lbStandardList, 'value');
        if any(iStd2Show),
            plotSpectra('StdFigure', iStd2Show);
        end
    end

%% UIcontrol callbacks: Reference deconvolution window
    function btnRefDeconvSpectrum_Callback (~,~)
        maxPPM=str2double(get(NMRdata.refDeconvData.txtRefDeconvMaxPPM, 'string'));
        minPPM=str2double(get(NMRdata.refDeconvData.txtRefDeconvMinPPM, 'string'));
        [~, minIndex]=min(abs(NMRdata.X.axisscale{2}-maxPPM));
        [~, maxIndex]=min(abs(NMRdata.X.axisscale{2}-minPPM));
        %Determine best peak based on minimum width and skewness
        [FWHH, skew]=getFWHH( NMRdata.X.data(:, (minIndex:maxIndex)), NMRdata.X.axisscale{2}(minIndex:maxIndex));
        qualityPeak=FWHH.*skew;
        [~, indexRef]=min(qualityPeak);
        set(NMRdata.refDeconvData.popRefDeconvSpectrum, 'value', indexRef)
    end

    function btnRefDeconv_Callback (~,~)
        if any(NMRdata.X.include{1}),
            NMRdata.XbackUp=NMRdata.X.data;
            acqPars=[];
            limits=zeros(1,2);
            limits(1)=str2double(get(NMRdata.refDeconvData.txtRefDeconvMaxPPM, 'string'));
            limits(2)=str2double(get(NMRdata.refDeconvData.txtRefDeconvMinPPM, 'string'));  
            aux=get(NMRdata.refDeconvData.popRefDeconvMode, 'value');
            listStr=get(NMRdata.refDeconvData.popRefDeconvMode, 'string');
            RefMode=listStr{aux};
            indexRef=get(NMRdata.refDeconvData.popRefDeconvSpectrum, 'value');
            switch RefMode
                case 'Signal in dataset'
                    RefMode=1;
                case 'Synthetic TMS signal'
                    RefMode=2;
            end
            indexes=getRang(NMRdata.X.axisscale{2}, limits(1), limits(2));
            [~, peaksPos]=max(real(NMRdata.X.data(:, indexes)), [], 2);
            newAxis=NMRdata.X.axisscale{2}(indexes);
            center=newAxis(mode(peaksPos));            
            [ NMRdata.X.data ] = referenceDeconv( NMRdata.X.data, NMRdata.X.axisscale{2}, center, limits, RefMode, NMRdata.X.SF(1), indexRef, NMRdata.X.include{1});
            close(NMRdata.hRefDeconvOptions);
            plotSpectra('MainFigure', NMRdata.visible);      
        end
    end

%% UIcontrol callbacks: reference shift
    function btnRefShif_Callback (~,~)
        if any(NMRdata.X.include{1}),
            NMRdata.XbackUp=NMRdata.X.data;
            tolerance=str2double(get(NMRdata.refShiftData.txtRefShifTolerance, 'string'));
            center=str2double(get(NMRdata.refShiftData.txtRefShifCenter, 'string'));
            NMRdata.X.data(NMRdata.X.include{1},:)=referenceShift(NMRdata.X.data(NMRdata.X.include{1},:), NMRdata.X.axisscale{2}, center, tolerance);
            close(NMRdata.hRefShiftOptions);
            plotSpectra('MainFigure', NMRdata.visible);  
        end
    end

%% UIcontrol callbacks: Align window
    function btnAlign_Callback(~,~)
        if any(NMRdata.X.include{1}),
            NMRdata.XbackUp=NMRdata.X.data;
            max_ppm=str2double(get(NMRdata.alignData.txtMaxAlignRegion, 'string'));
            min_ppm=str2double(get(NMRdata.alignData.txtMinAlignRegion, 'string'));
            type=get(NMRdata.alignData.cbApplyAllSpectrum, 'value');
            if type==1,
                type='all';
            else
                type='range';
            end
            NMRdata.X.data(NMRdata.X.include{1},:)=alignbyRegion(NMRdata.X.data(NMRdata.X.include{1},:), NMRdata.X.axisscale{2}, max_ppm, min_ppm, type);
            close(NMRdata.hAlignOptions);
            plotSpectra('MainFigure', NMRdata.visible);      
        end
    end

%% UIcontrol callbacks: Phasing window
    function btnAddRange_Callback(hObject,~)
        hParentName=get(hObject, 'tag');
        oldRegions = get(NMRdata.(hParentName).tblRanges,'Data');
        newRegions=zeros(1, 2);
        newRegions = [oldRegions; newRegions];
        set(NMRdata.(hParentName).tblRanges, 'Data',newRegions)
    end

    function btnDelRange_Callback(hObject,~)
        hParentName=get(hObject, 'tag');
        if isfield(NMRdata.(hParentName), 'delRow'),
            oldRegions = get(NMRdata.(hParentName).tblRanges,'Data');

            if size(oldRegions)>1,
                oldRegions(NMRdata.(hParentName).delRow,:)=[];
            end
            set(NMRdata.(hParentName).tblRanges,'Data',oldRegions)
        end
    end

    function tblRanges_CellSelectionCallback(hObject,eventdata)
        hParentName=get(hObject, 'tag');
        NMRdata.(hParentName).delRow=eventdata.Indices(:,1);
    end

    function btnPhase_Callback(~,~)
        if any(NMRdata.X.include{1}),
            NMRdata.XbackUp=NMRdata.X.data;
            mode=get(NMRdata.phaseData.popPhaseMode, 'value');
            set(NMRdata.phaseData.btnPhase,  'Enable', 'off');
            NMRdata.phaseRegions=get(NMRdata.phaseData.tblRanges,'Data');
            auto=get(NMRdata.phaseData.cbPhaseAutoRegions, 'value');
            [NMRdata.X.data(NMRdata.X.include{1},:)] = phase1D( NMRdata.X.data(NMRdata.X.include{1},:), NMRdata.X.axisscale{2}, mode, NMRdata.phaseRegions, 0.3, 10, auto);
            [NMRdata.X.data(NMRdata.X.include{1},:)] = phase1D( NMRdata.X.data(NMRdata.X.include{1},:), NMRdata.X.axisscale{2}, 0, [], 0.3, [], 1);
            close(NMRdata.hPhaseOptions);
            plotSpectra('MainFigure', NMRdata.visible);
        else
          warndlg('No sample(s) selected');
        end
    end

    function cbAutoRange_Callback(hObject,~)
        hParentName=get(hObject, 'tag');
        aux=get(hObject, 'value');
        if aux==1,
            set(NMRdata.(hParentName).tblRanges, 'Data',zeros(1, 2))
            set(NMRdata.(hParentName).tblRanges, 'Enable', 'off');
            set(NMRdata.(hParentName).btnAddRange, 'Enable', 'off');
            set(NMRdata.(hParentName).btnDelRange, 'Enable', 'off');
        else
            set(NMRdata.(hParentName).tblRanges, 'Enable', 'on');
            set(NMRdata.(hParentName).btnAddRange, 'Enable', 'on');
            set(NMRdata.(hParentName).btnDelRange, 'Enable', 'on');
        end         
    end

    function popPhaseMode_Callback (hObject,~)
        mode=get(hObject, 'value');
        if mode==2,
            set(NMRdata.phaseData.cbPhaseAutoRegions, 'Enable', 'on');
        else
            set(NMRdata.phaseData.cbPhaseAutoRegions, 'Enable', 'off');
            set(NMRdata.phaseData.cbPhaseAutoRegions, 'value', 0);
            set(NMRdata.phaseData.tblRanges, 'Enable', 'on');
            set(NMRdata.phaseData.btnAddRange, 'Enable', 'on');
            set(NMRdata.phaseData.btnDelRange, 'Enable', 'on');            
        end
    end

%% UIcontrol callbacks: Baseline window
    function popBaselineMode_Callback(hObject, ~)
        str=get(hObject, 'string');
        aux=get(hObject, 'value');
        switch str{aux}
            case {'median substraction','cubic spline', 'cubic Hermite'}
                set (NMRdata.baselineData.lblOrder, 'visible', 'off');
                set (NMRdata.baselineData.txtOrder, 'visible', 'off');                
            case 'polynomial'
                set (NMRdata.baselineData.lblOrder, 'visible', 'on');
                set (NMRdata.baselineData.txtOrder, 'visible', 'on');
        end
    end

    function btnBaseline_Callback(~,~)
        if any(NMRdata.X.include{1}),
            NMRdata.XbackUp=NMRdata.X.data;
            str=get(NMRdata.baselineData.popBaselineMode, 'string');
            aux=get(NMRdata.baselineData.popBaselineMode, 'value');   
            binSize=0.2;
            switch str{aux}
                case 'median substraction'
                    binSize=0.8;
                    mode=0;
                case 'cubic spline'
                    mode=1;
                case 'cubic Hermite'
                    mode=2;
                case 'polynomial'
                    mode=3;
            end 
            order=str2double(get(NMRdata.baselineData.txtOrder, 'string'));
            set(NMRdata.baselineData.btnBaseline,  'Enable', 'off');
            NMRdata.bcRegions=get(NMRdata.baselineData.tblRanges,'Data');
            regions=NMRdata.bcRegions;
            if get(NMRdata.baselineData.cbBaselineAutoRanges, 'value'),
                regions=[];
            end
            [NMRdata.X.data(NMRdata.X.include{1},:)] = baselineCorrection(NMRdata.X.data(NMRdata.X.include{1},:), NMRdata.X.axisscale{2}, mode, regions, binSize, order);
            close(NMRdata.hBaselineOptions);
            plotSpectra('MainFigure', NMRdata.visible);     
        end
    end

%% UIcontrol callbacks: Import window
    function btnImport_Callback(~,~)
        %Load data from uicontrol
        NMRdata.importData.options_import.files_path=get(NMRdata.importData.txtPath, 'string');
        NMRdata.importData.options_import.file_filter=get(NMRdata.importData.txtFilter, 'string');         
        NMRdata.importData.options_import.exp_no=str2double(get(NMRdata.importData.txtExp_no, 'string'));
        if isfield(NMRdata.importData, 'txtProc_no'),
            NMRdata.importData.options_import.proc_no=str2double(get(NMRdata.importData.txtProc_no, 'string'));
        end
        %Add additional data needed
        NMRdata.importData.options_import.minppm=-50;
        NMRdata.importData.options_import.maxppm=50;
        %Launch import routine
        [NMRdata.X, exitFlag]=importBrukerBatch(NMRdata.importData.options_import);
        if exitFlag > 0
            return;
        end
        NMRdata.X.P1=NMRdata.X.P1*1e-6;
        NMRdata.X.SF=NMRdata.X.SF*1e6;
        if NMRdata.importData.typeFile==2,
            %Apply FT
            NMRdata.X.include{1}=true(size(NMRdata.X.data,1),1);
            %First shift the spectra according to bruker delay GRPDLY
            %and set last points to zero
            NMRdata.X.data=circshift(NMRdata.X.data, [0 -round(NMRdata.X.GRPDLY)]);
            NMRdata.X.data(:, end-round(NMRdata.X.GRPDLY):end)=0;
            %Apply windowing
            if ~strcmp(NMRdata.importData.WDWtype, 'non');
                lb=str2double(get(NMRdata.importData.txtWDWF2lb, 'string'));
                gb=str2double(get(NMRdata.importData.txtWDWF2gb, 'string'));
                f=getApodizationFunction(NMRdata.importData.WDWtype, NMRdata.X.axisscale{2}, lb, gb, 0);
            else
                f=ones(1,length(NMRdata.X.axisscale{2}));
            end
            %Apply offset reduction
            if get(NMRdata.importData.cbFirstPointby05, 'value'),
                NMRdata.X.data(:,1)=NMRdata.X.data(:,1)*0.5;
            end
            f=repmat(f, size(NMRdata.X.data,1), 1);
            if size(NMRdata.X.data,2)>size(f,2),
                %In case file size differ from TD
                NMRdata.X.data(size(f,2)+1:end)=[];
            end                
            NMRdata.X.data=NMRdata.X.data.*f;
            NMRdata.X.SI=NMRdata.X.TD/2;
            zeroFillingF2=str2double(get(NMRdata.importData.txtZerofillingF2, 'string'));
            if zeroFillingF2>0,
                value=nextpow2(NMRdata.X.SI(1));
                NMRdata.X.SI=2^(value+zeroFillingF2);
                NMRdata.X.data=[NMRdata.X.data zeros(size(NMRdata.X.data,1), NMRdata.X.SI-size(NMRdata.X.data,2))];
            end
            NMRdata.X.data=fliplr(fftshift(fft(NMRdata.X.data', NMRdata.X.SI(1)),1)');
            NMRdata.X.data=real(NMRdata.X.data)+1*sqrt(-1)*imag(NMRdata.X.data);
            offsetPPM=NMRdata.X.O1(1)*1e6/NMRdata.X.SF(1);
            NMRdata.X.OFFSET(1)=offsetPPM;
            SW1PPM=NMRdata.X.SW1(1)*1e6/NMRdata.X.SF(1);
            NMRdata.X.axisscale{2}=linspace((SW1PPM/2)+offsetPPM, -(SW1PPM/2)+offsetPPM, NMRdata.X.SI(1)); 
        end
        NMRdata.scaleXaxis{2}=NMRdata.X.axisscale{2};
        NMRdata.X.include{1}=1:size(NMRdata.X.data,1);
        NMRdata.visible=NMRdata.X.include{1};
        included=num2cell(true(size(NMRdata.X.data,1),1));
        cellArray=[included, NMRdata.X.label{1}];
        set(NMRdata.uiControls.tblIncludedSamples, 'data', cellArray);
        set(NMRdata.uiControls.tblIncludedSamples, 'Units', 'pixels');
        units=max(max(cellfun('size', NMRdata.X.label{1}, 2), 6));
        fontsize=get(NMRdata.uiControls.tblIncludedSamples, 'fontsize');
        set(NMRdata.uiControls.tblIncludedSamples,'ColumnWidth', {'auto'  units*fontsize});
        set(NMRdata.uiControls.lblParResult, 'string', '');        
        for j=1:6,
            set(NMRdata.uiControls.tblresults(j), 'Rowname', NMRdata.X.label{1});
            set(NMRdata.uiControls.tblresults(j), 'columNname', []);
            set(NMRdata.uiControls.popCorrectWidth, 'String', {'None'});
            set(NMRdata.uiControls.popCorrectWidth, 'value', 1);
            set(NMRdata.uiControls.tblresults(j), 'data', []);
        end
        if NMRdata.fitPars.normAcqPars==1,
            NMRdata.X.data=NMRdata.X.data./repmat((NMRdata.X.RG'.*NMRdata.X.NS'.*NMRdata.X.P1'), [1 size(NMRdata.X.data, 2) size(NMRdata.X.data, 3)]);
        end
        if isfield(NMRdata, 'standards'),
            NMRdata=rmfield(NMRdata, 'standards');
        end  
        LoadFileStandard(NMRdata.paths.stds, []);
        set(NMRdata.uiControls.btnSelAll, 'enable', 'on');
        set(NMRdata.uiControls.lbVisible, 'enable', 'on');
        set(NMRdata.uiControls.lbVisible, 'value', []);
        set(NMRdata.uiControls.lbVisible, 'string', [{'All'}; NMRdata.X.label{1}]);
        NMRdata.XbackUp=[];
        close(NMRdata.hImportOptions);
        plotSpectra('MainFigure', NMRdata.X.include{1});        
    end

    function searchDirPath_Callback(~, ~)
        filePath = uigetdir;
        set(NMRdata.importData.txtPath, 'string', filePath);
    end

    function popWDW_Callback(hObject, ~)
        aux=get(hObject, 'value');
        listStr=get(hObject, 'string');
        NMRdata.importData.WDWtype=listStr{aux};
        NMRdata.importData.WDWtype(4:end)=[];
        switch NMRdata.importData.WDWtype
            case 'non'
                set(NMRdata.importData.txtWDWF2lb, 'string', '');
                set(NMRdata.importData.txtWDWF2lb, 'enable', 'off');
                set(NMRdata.importData.txtWDWF2gb, 'string', '');
                set(NMRdata.importData.txtWDWF2gb, 'enable', 'off');                
            case 'lor'
                set(NMRdata.importData.txtWDWF2lb, 'enable', 'on');
                set(NMRdata.importData.txtWDWF2lb, 'string', '0');
                set(NMRdata.importData.txtWDWF2gb, 'string', '');
                set(NMRdata.importData.txtWDWF2gb, 'enable', 'off');                  
            case 'gau'
                set(NMRdata.importData.txtWDWF2lb, 'enable', 'on');
                set(NMRdata.importData.txtWDWF2lb, 'string', '0');
                set(NMRdata.importData.txtWDWF2gb, 'enable', 'on');
                set(NMRdata.importData.txtWDWF2gb, 'string', '0');                
        end
    end

    function mnPath_Callback(~, ~, sourceHandle)
        cAux=clipboard('paste');
        if ischar(cAux) && (strcmp(cAux(2:3), ':\') || strcmp(cAux(1), '/')),
           set(sourceHandle, 'string', cAux);
        end
    end
    
%% Auxiliary functions   
    function QuitLipspin(~,~)
        if isfield(NMRdata,'hImportOptions')
            close(NMRdata.hImportOptions)
        end
        if isfield(NMRdata,'hPhaseOptions')
            close(NMRdata.hPhaseOptions)
        end        
        delete(hMainFigure)
    end

    function SavePlot(hObject,~)
        [fileName, filepath] = uiputfile('*.png', 'Save file', '');
        h1=get(hObject, 'Parent');
        h2=get(h1, 'Parent');
        tag=get(h2, 'tag');
        switch tag
            case 'StdFigure'
                h=NMRdata.loadStandardsData.hStandardAxes;
            case 'MainFigure'
                h=NMRdata.hMainAxes;
        end
        if ~isnumeric(fileName) && isfield(NMRdata, 'hPlot'),
            hAux=figure('Visible', 'off');
            copyobj(h, hAux);
            set(gca, 'Units','normalized', 'Position',[0.05 0.05 0.9 0.9]);
            xlbl=get(NMRdata.lblMainXaxis, 'string');
            xlabel(xlbl);
            print(hAux, [filepath fileName], '-dpng', '-r300');
            close(hAux);
        end
    end

    function plotSpectra(type, visible)
        zoom off
        pan off
        datacursormode off
        switch type
            case 'MainFigure'
                ylimits=[min(min(real(NMRdata.X.data(visible,:)))) max(max(real(NMRdata.X.data(visible,:))))];
                NMRdata.hPlot=plot(NMRdata.hMainAxes, NMRdata.scaleXaxis{2}, real(NMRdata.X.data(visible,:)));
                linelabel(NMRdata.hPlot, NMRdata.X.label{1}(visible,:));
                ylim(ylimits);
        %             axis tight
                zoom reset
                xlim([NMRdata.scaleXaxis{2}(end) NMRdata.scaleXaxis{2}(1)]);
                set(NMRdata.hMainAxes, 'xdir', 'reverse');                
            case 'StdFigure'
                if isempty(NMRdata.standards.data),
                    legend(NMRdata.loadStandardsData.hStandardAxes, 'hide','Interpreter', 'none')
                    cla(NMRdata.loadStandardsData.hStandardAxes);
                else
                    plot(NMRdata.loadStandardsData.hStandardAxes, NMRdata.X.axisscale{2}, NMRdata.standards.data(visible,:));
                    set(NMRdata.loadStandardsData.hStandardAxes, 'xdir', 'reverse');   
                    pos=get(NMRdata.loadStandardsData.hStandardAxes, 'position');
                    pos(1)=pos(1)+pos(3)*0.85;
                    pos(3)=pos(3)*0.10;
                    h=legend(NMRdata.loadStandardsData.hStandardAxes, NMRdata.standards.label{1}(:,visible), 'location', 'east', 'position', pos, 'fontsize', 7 , 'Interpreter', 'none');
                    set(h,'PlotBoxAspectRatioMode','manual');
                    set(h,'PlotBoxAspectRatio',[1 0.8 1]); 
                end
        end
        datacursormode on
    end

    function personaliseToolbar(h)
        tag=get(h, 'tag');
        hToolbar = findall(h,'tag','FigureToolBar');
        hTools=allchild(hToolbar);
        set(hTools,'Visible','off');
        set(hTools,'Separator','off');
        
        path=which('lipspin');
        
        cdataSaveFigure = imread('Images/SaveFigure.png');
        uipushtool(hToolbar, 'cdata',cdataSaveFigure, 'tooltip','Save Figure', 'ClickedCallback',{@SavePlot},'Visible','on','Separator','on');
        if strcmp(tag, 'MainFigure'),
            cdataUndoX = imread('Images/Undo.png');
            uipushtool(hToolbar, 'cdata', cdataUndoX, 'tooltip','Undo last action', 'ClickedCallback',{@UndoX},'Visible','on','Separator','on');
        end
        
        aux = findall(hTools,'tag','Exploration.ZoomIn');
        set(aux,'Visible','on');
        aux = findall(hTools,'tag','Exploration.ZoomOut');
        set(aux,'Visible','on');        
        aux = findall(hTools,'tag','Exploration.Pan');
        set(aux,'Visible','on');
        aux = findall(hTools,'tag','Exploration.DataCursor');
        set(aux,'Visible','on');
    end

    function LoadFileStandard(pathname, filename)
        if exist(NMRdata.paths.stds, 'file')
            if isempty(filename),
                stds=dir([NMRdata.paths.stds '\*.nmrstd']);
                pathname=[NMRdata.paths.stds '\'];
            else
                stds.name=filename;
            end
            if ~isempty(stds),
                nbStds=size(stds,1);
                for i=1:nbStds,
                    f=fopen([pathname '\' stds(i).name]);
                    patternFile = fread(f, '*char')';
                    fclose(f);
                    if ~strcmp(patternFile,'')
                        line=1;
                        posCR=[0 strfind(patternFile, sprintf('\n'))];
                        while (line<length(posCR))
                            currentLine=patternFile(posCR(line)+1:posCR(line+1));
                            aux=strfind(currentLine, '=');
                            switch currentLine(3:aux(1)-1)          
                                case 'NAME'
                                    name=(currentLine(aux+1:end-1));
                                case 'LOWLIMIT'
                                    lowlimit=str2double(currentLine(aux+1:end-1));
                                case 'HIGHLIMIT'
                                    highlimit=str2double(currentLine(aux+1:end-1));
                                case 'FREQUENCY'
                                    sf1=str2double(currentLine(aux+1:end-1));
                                case 'SPECTRUM'
                                    f = fopen([pathname stds(i).name]);
                                    fseek(f, posCR(line+1), -1);
                                    spec = fread(f,'double')';
                                    fclose(f);
                                    line=length(posCR);
                            end
                            line=line+1;
                        end
                        axisscale=highlimit:-(highlimit-lowlimit)/(length(spec)-1):lowlimit;
                    end
                    if axisscale(1)<NMRdata.X.axisscale{2}(1)
                        axisscale=[NMRdata.X.axisscale{2}(1) axisscale];
                        spec=[0 spec];
                    else
                        [~, ind]=min(abs(axisscale-NMRdata.X.axisscale{2}(1)));
                        if axisscale(ind)>NMRdata.X.axisscale{2}(1),
                            ind=ind+1;
                        end
                        axisscale(1:ind)=[];
                        spec(1:ind)=[];
                    end
                    if axisscale(end)>NMRdata.X.axisscale{2}(end)
                        axisscale=[axisscale NMRdata.X.axisscale{2}(end)];
                        spec=[spec 0];
                    else
                        [~, ind]=min(abs(axisscale-NMRdata.X.axisscale{2}(end)));
                        if axisscale(ind)<NMRdata.X.axisscale{2}(end),
                            ind=ind-1;
                        end
                        axisscale(ind:end)=[];
                        spec(ind:end)=[];
                    end 
                    spec = interp1(axisscale, spec, NMRdata.X.axisscale{2});
                    newInd=1;
                    if isfield(NMRdata, 'standards')
                        newInd=size(NMRdata.standards.data, 1)+1;
                    end
                    spec(isnan(spec))=0;
                    NMRdata.standards.data(newInd,:)=spec;
                    NMRdata.standards.label{1}{newInd}=name;
                end
            [NMRdata.standards.label{1} index]=sort(NMRdata.standards.label{1});
            NMRdata.standards.data=NMRdata.standards.data(index,:);
            else
                warndlg('Standards not found');
            end
        else
            warndlg('Path for standard spectra not found');
        end
    end

    function EnableDisableUIcontrols(handlers, state)
        for i=1:length(handlers)
            set(handlers(i), 'Enable', state);
        end
    end

    function UndoX(~,~)
        NMRdata.X.data=NMRdata.XbackUp;
        plotSpectra('MainFigure', NMRdata.visible);
    end
end






