function [data] =readBrukerFile(filename, fileType, byteorder, dataDim, memoryDim)
%     Read bruker raw 1D,2D FID and 2D FT spectra (tested with Topspin 2.x & 3.x versions and bruker Avance spectrometers)
%
%     INPUT
%     filename: full file path
%     filetype: 
%              'FID' - time domain decaying signal 
%              'SPEC' - FT spectrum
%     byteorder: 
%              'l' - little endian data storage format
%              'b' - big endian data storage format
%     dataDim: column vector containing dimensions of F2, F1 axes (only for 2D data, otherwise [])
%     memoryDim: storage data blocks in memory (only for 2D data, see bruker manuals for details)
%
%     OUTPUT
%     data: raw NMR data

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

dataformat = 'l';
if byteorder==1
    dataformat = 'b';
end;

dim=[];

switch fileType
    case 'FID'
        if exist([filename 'fid'], 'file')
            filename = [filename 'fid'];
            dim=1;
        elseif exist([filename 'ser'], 'file')
            filename = [filename 'ser'];
            dim=2;
        else
            fprintf('File not found: %s\n \r', filename );
            return;
        end
        %Load file
        [f, MESSAGE] = fopen(filename,'r',dataformat);
        data = fread(f,'int32');
        fclose(f);
    case 'SPEC'
        if exist([filename '2rr'], 'file') && exist([filename '2ii'], 'file'),
            realfile = [filename '2rr'];
            imagfile = [filename '2ii'];
            dim=2;
        elseif exist([filename '1r'], 'file') && exist([filename '1i'], 'file'),
            realfile = [filename '1r'];
            imagfile = [filename '1i'];
            dim=1;
        elseif exist([filename '2rr'], 'file'),
            fprintf('Warning: imag file not found');
            realfile = [filename '2rr'];           
            dim=2;
        elseif exist([filename '1r'], 'file'),
            fprintf('Warning: imag file not found');
            realfile = [filename '1r'];
            dim=1;
        else    
            fprintf('File not found: %s\n \r', filename );
            return;
        end
        %Load file
        [f, MESSAGE] = fopen(realfile,'r',dataformat);
        real = fread(f,'int32');
        fclose(f);
        imag=[];
        if exist('imagfile'),
            [f, MESSAGE] = fopen(imagfile,'r',dataformat);
            imag = fread(f,'int32');
            fclose(f);
        end
end

%rearrange data depending on the type of file loaded
if strcmp(fileType, 'FID'),
    data=data(1:2:end)-1*sqrt(-1)*data(2:2:end);
    if dim==2,
        dataDim1=dataDim/2;
        dataDim2=length(data)/dataDim1;
        data=reshape(data, dataDim1, dataDim2);
    end
else
    if ~isempty(imag),
        data=real+sqrt(-1)*imag;
    else
        data=real;
    end
    if dim==2,
        %if 2D data, re-arrangement of spectrum file is neccessary
        NoSM = dataDim(1)*dataDim(2)/(memoryDim(1)*memoryDim(2));    % Total number of Submatrixes
        NoSM2 = dataDim(2)/memoryDim(2);               % No of SM along F1    
        data=reshape(permute(reshape(permute(reshape(real,memoryDim(1),memoryDim(2),NoSM), [2 1 3]), memoryDim(2),dataDim(1),NoSM2), [2 1 3]),dataDim(1),dataDim(2))';
    end
end

end



