function [out, sample_names, exitFlag] = readBrukerParameter (options, param, subparam, listfolders)
%   Read one bruker parameter for multiple 1D or 2D NMR samples at a time
% 
%   INPUT:
%   options: file(s) options
%            files_path: Sample location
%            file_filter: List only samples containing the filter
%            exp_no: experimental folder name (as a natural number)
%            type: 'acqus' (FID files) or 'procs' (FT spectra files)
%   param: parameter ID
%   subparam: number for vector parameters (In case a number follows the parameter, e.g. D1-D20)
%   listfolders: [1xN] cell of strings with sample names (allows skipping search window with listdlg)
%
%   OUTPUT:
%   out: [1xN] row vector with parameter values. N=number of samples
%   sample_names: names for samples with found parameter files
%   exitFlag: exception flag

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

out=[];

% Check coherence
exitFlag=0;
if ~exist(options.files_path, 'dir'), 
    exitFlag=1;
    disp('Directory does not exist');
    return;   
end

if nargin==4,
    sample_folders=listfolders;
else
    d=dir(options.files_path);
    d = d(find(cellfun(@(x)(x==1),{d(:).isdir})));
    %susbtract '.' and '..' folders
    d = d(find(~cellfun(@(x)(strncmp('.', x, 1)), {d(:).name})));
    if ~isempty(options.file_filter)
        d = d(find(~cellfun(@(x)(isempty((strfind(x, options.file_filter)))), {d(:).name})));
    end
    folder_list = {d.name};
    [index, ok] = listdlg('ListString', folder_list);
    if ok == 0
        return;
    end
    sample_folders = folder_list(index);
end

sample_names = sample_folders;
i=1;
while i<=length(sample_names)
    if strcmp(options.files_path(end), '\'),
        options.files_path(end)=[];
    end    
    if strcmp(options.type, 'acqus') || strcmp(options.type, 'acqu2s'), 
        sample_folders{i} = [options.files_path '\' sample_folders{i} '\' num2str(options.exp_no) '\'];
    else
        sample_folders{i} = [options.files_path '\' sample_folders{i} '\' num2str(options.exp_no) '\pdata\' num2str(options.proc_no) '\'];
    end
    if ~isdir(sample_folders{i})
        sample_names(i)=[];
        sample_folders(i)=[];
    else
        i=i+1;
    end
end    

aux2=subparam;    

for i=1:length(sample_folders),
    subparam=aux2;
    fn = [sample_folders{i} options.type];
    TXT = LoadCompleteText(fn);
    if (strcmp(TXT,'') == 1)
        exitFlag=2;
        disp('Parameters file is void');
        return;   
    end;
    POS=[0 strfind(TXT,sprintf('\n'))];    
    fpos=1;
    while (fpos<length(POS))
       [pp, fpos] = GetTextLine(TXT,POS,fpos);
       if (strncmp(pp,['##$' param '= '],5+length(param)) == 1) 
           if ~isempty(subparam),
               [pp, fpos] = GetTextLine(TXT,POS,fpos);
               aux=find(pp==' ');
               while (subparam > length(aux)),
                   subparam=subparam-length(aux)-1;
                   [pp, fpos] = GetTextLine(TXT,POS,fpos);
                   aux=find(pp==' ');
               end
               if (subparam ==0),
                   out(i) = str2double(pp(1:aux(subparam+1)));
               else
                   out(i) = str2double(pp(aux(subparam):aux(subparam+1)));
               end
           else
               out(i) = str2double(pp(6+length(param):length(pp)-1));
           end
       end
    end
end
if isempty(out),
    exitFlag=3;
    disp(['Parameter ' param ' not found']);
    return;       
end


function [Line,fposnew] = GetTextLine(TXT,POS,fpos)
    Line=TXT(POS(fpos)+1:POS(fpos+1));
    fposnew=fpos+1;

function [TXT, TXTO] = LoadCompleteText(filename)
    f = fopen(filename); 
    if f<0 
       fprintf('Error: file "%s" not found\n', filename);
       TXT = '';
       return;
    end;
    TXTO = fread(f);
    fclose(f);
    TXT = char(TXTO');