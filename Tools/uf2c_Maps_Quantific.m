function varargout = uf2c_Maps_Quantific(varargin)
% UF²C M-file for uf2c_Maps_Quantific.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2017
%
% Copyright (c) 2017, Brunno Machado de Campos
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_Maps_Quantific_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_Maps_Quantific_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function uf2c_Maps_Quantific_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_Maps_Quantific_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddMaps_Callback(hObject, eventdata, handles)
global file path
[file,path] = uigetfile('*.nii','Add the subjects MAPS','MultiSelect','on');
if ~iscell(file)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    file = {file};
end
set(handles.text2,'String',[num2str(numel(file)) ' Maps added'])

function AddROIs_Callback(hObject, eventdata, handles)
global fileR pathR

if get(handles.AddR,'Value')
    [fileR,pathR] = uigetfile('*.nii','Add all ROIs Masks','MultiSelect','on');
    if ~iscell(fileR)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileR= {fileR};
    end
    set(handles.text3,'String',[num2str(numel(fileR)) ' ROIs added'])
else
    [fileR,pathR] = uigetfile('*.nii','Add the Mask','MultiSelect','off');
    if ~iscell(fileR)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileR = {fileR};
    end
    set(handles.text3,'String','Mask added')
end

function AddR_Callback(hObject, eventdata, handles)
if get(handles.AddR,'Value')
    set(handles.AddM,'Value',0)
    set(handles.AddROIs,'String','Add all ROIs')
else
    set(handles.AddM,'Value',1)
    set(handles.AddROIs,'String','Add a mask')
end

function AddM_Callback(hObject, eventdata, handles)
if get(handles.AddM,'Value')
    set(handles.AddR,'Value',0)
    set(handles.AddROIs,'String','Add a mask')
else
    set(handles.AddR,'Value',1)
    set(handles.AddROIs,'String','Add all ROIs')
end

function RunB_Callback(hObject, eventdata, handles)
global fileR pathR file path OutDir

set(handles.text6,'String','Running...')
drawnow
fprintf('Process started === %s\r\n',spm('time'));

% dateNow = clock;
% Fname = (['1-RESULTs_MapQuant_' num2str(dateNow(3)) '_' num2str(dateNow(2)) '_' num2str(dateNow(1)) '--' num2str(dateNow(4)) '_' num2str(dateNow(5)) '.txt']);
verifMAP = nifti([path file{1}]);

if ~get(handles.wiq,'Value')
    verifROI = nifti([pathR fileR{1}]);
    Thre = get(handles.MaskThre,'String');
    ClusThre = str2num(get(handles.edit4,'String'));
    if ~isequal(verifROI.dat.dim(1:3),verifMAP.dat.dim(1:3))
        warndlg({'Your ROIs and your Maps have different dimenions. UF²C requires';...
            'files with same dimension to guarantee high quality results.'},'Attention!');
        return
    end
else
    ThresholdMap = 0;
end

ThresholdMap = str2num(get(handles.MapThre,'String'));

if isempty(OutDir)
    OutDir = [path filesep];
else
    OutDir = [OutDir filesep];
end

fideGx1 = fopen([OutDir,get(handles.Resname,'String'),'.txt'],'w+'); % CHARGE OUTPUT LOG FILE

if ~get(handles.wiq,'Value')
    if get(handles.AddM,'Value')
        a = nifti([pathR fileR{1}]);
        matx = a.dat(:,:,:,1);
        if ~isequal(Thre,'none') % if is not binarized, the code will binarize and runs
            if isequal(numel(unique(matx)),2) % sometimes.... a "not" binarized image IT IS binarized...
                L = bwlabeln(matx);

                vet = unique(L);
                vet = vet(2:end);

                iidx = 0;
                for i = 1:size(vet,1)
                    eval(sprintf('Mat%d = double((L==%d));',i,vet(i)));
                    if nnz(eval(sprintf('Mat%d',i)))>ClusThre
                        iidx = iidx+1;
                        struQ = a;
                        struQ.dat.fname = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                        struQ.dat.dtype = 'FLOAT32-LE';
                        eval(sprintf('struQ.dat(:,:,:) = Mat%d;',i));
                        create(struQ)
                        clear struQ
                        filesR{iidx,1} = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                    end
                end
            else
                matx = matx>str2num(Thre);
                L = bwlabeln(matx);
                vet = unique(L);
                vet = vet(2:end);
                iidx = 0;
                for i = 1:size(vet,1)
                    eval(sprintf('Mat%d = double((L==%d));',i,vet(i)));
                    if nnz(eval(sprintf('Mat%d',i)))>ClusThre
                        iidx = iidx+1;
                        struQ = a;
                        struQ.dat.fname = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                        eval(sprintf('struQ.dat(:,:,:) = Mat%d;',i));
                        create(struQ)
                        filesR{iidx,1} = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                    end
                end
            end
        else
            if ~isequal(numel(unique(matx)),2) % Sometimes.... a "binarized" image IT IS NOT......
                fprintf('Image %s is not binarized...\r\n',fileR{1})
                fprintf('Defines a Threshold in a next time!\r\n')
                beep
                beep
                return
            else
                L = bwlabeln(matx);
                vet = unique(L);
                vet = vet(2:end);
                iidx = 0;
                for i = 1:size(vet,1)
                    eval(sprintf('Mat%d = double((L==%d));',i,vet(i)));
                    if nnz(eval(sprintf('Mat%d',i)))>ClusThre
                        iidx = iidx+1;
                        struQ = a;
                        struQ.dat.fname = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                        eval(sprintf('struQ.dat(:,:,:) = Mat%d;',i));
                        create(struQ)
                        filesR{i,1} = sprintf('%s%d-%s',OutDir,iidx,fileR{1});
                    end
                end
            end
        end
        fileR = filesR;
        fileR = fileR(~cellfun(@isempty, fileR));
        set(handles.text3,'String',[num2str(numel(filesR)) ' splitted clusters'])

    end

    for il = 1:numel(fileR)
        if  get(handles.AddR,'Value')
            fileR{il} = [pathR,fileR{il}];
        else
            [~,bx,cx] = fileparts(fileR{il});
            fileR{il} = [OutDir,bx,cx];
        end
    end
    
    try
        [roi12] = uf2c_AnatDescrip(fileR,OutDir,[get(handles.Resname,'String'),'_AnatDescrip.txt'],'mask');
    catch
        for il = 1:numel(fileR)
            tmpNS = nifti(fileR{il});
            roi12(:,:,:,il) = tmpNS.dat(:,:,:);
        end
    end

    result = cell(size(file,2),size(fileR,2)+1);
    result{1,1} = 'Name';
    fprintf(fideGx1,'Name \t');

    for i2 = 2:numel(fileR)+1
        [~,bx,~] = fileparts(fileR{i2-1});
        result{1,i2} = bx;
        fprintf(fideGx1,'%s \t ',bx);
    end
end
fprintf(fideGx1,'Whole mask/ALL ROIs AVG');
fprintf(fideGx1,'\r\n');

for i = 1:numel(file)
    file2{i} = [path file{i}];
end
verifMAP = nifti(file2{1});
if isequal(numel(verifMAP.dat.dim),4)
    warndlg('One or more maps added has a 4th dimension. The first volume of the 4th dimension will be considered.','Attention!');
end

for i = 1:numel(file2)
   [~,a2,~] = fileparts(file2{i});
   result{i+1,1} = a2;
   fprintf(fideGx1,'%s \t ',a2);
   
   stru = nifti(file2{i});
   matI = double(stru.dat(:,:,:,1));
   matI(matI(:,:,:,1)<ThresholdMap) = 0;

   if get(handles.wiq,'Value')
       rela1 = sum(sum(sum(matI)));
       nofVoxs = nnz(matI);
       rela1 = rela1/nofVoxs;
       
       fprintf(fideGx1,'%.4f\t  NofVoxels:\t%d',rela1,nofVoxs);
       fprintf(fideGx1,'\r\n');
   else
        for j = 1:numel(fileR)
            mat1 = double(roi12(:,:,:,j));
            if ~isequal(Thre,'none')
                mat1(mat1<str2num(Thre)) = 0;
            end
            try
                ResMat1 = matI.*mat1;
                ResMat1(isnan(ResMat1)) = 0;
                catch
                    warndlg('Maybe, the  matrices sizes of the Maps and ROIs do not matches','Process aborted!')
                    return
            end
            nume1 = sum(sum(sum(ResMat1)));
            nz1 = nnz(ResMat1);
            rela1 = nume1/nz1;
            result{i+1,j+1} = rela1;
            fprintf(fideGx1,'%.4f\t ',rela1);
            clear mat
        end

        clear matl
        %To add an average of all clusters/ROIs
        mat1 = sum(roi12,4);
        mat1(mat1>1) = 1;

        if ~isequal(Thre,'none')
            mat1(mat1<str2num(Thre)) = 0;
        end
        try
            ResMat1 = matI.*mat1;
            ResMat1(isnan(ResMat1)) = 0;
            catch
                warndlg('Maybe, the  matrices sizes of the Maps and ROIs do not matches','Process aborted!')
                return
        end
        nume1 = sum(sum(sum(ResMat1)));
        nz1 = nnz(ResMat1);
        rela1 = nume1/nz1;
        result{i+1,j+2} = rela1;
        fprintf(fideGx1,'%.4f\t ',rela1);
        %%%%%

        fprintf(fideGx1,'\r\n');
   end
end

fclose('all');
set(handles.text6,'String','Done!')
fprintf('Done! ============== %s\r\n',spm('time'));

function OutDir_Callback(hObject, eventdata, handles)
global OutDir

OutDir = uigetdir('','Define the output directory');
if isequal(OutDir,0)
    clear OutDir;
else
    set(handles.text9,'String',OutDir)
end

function Resname_Callback(hObject, eventdata, handles)

function MapThre_Callback(hObject, eventdata, handles)

function MaskThre_Callback(hObject, eventdata, handles)

function MaskThre_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MapThre_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function wiq_Callback(hObject, eventdata, handles)
if get(handles.wiq,'Value')
    set(handles.AddM,'Value',0)
    set(handles.AddR,'Value',0)
    set(handles.AddM,'Enable','off')
    set(handles.AddR,'Enable','off')
    set(handles.text10,'Enable','off')
    set(handles.text4,'Enable','off')
    set(handles.edit4,'Enable','off')
    set(handles.MaskThre,'Enable','off')
    set(handles.AddROIs,'Enable','off')
else
    set(handles.AddM,'Value',0)
    set(handles.AddR,'Value',1)
    set(handles.AddROIs,'Enable','on')
    set(handles.AddM,'Enable','on')
    set(handles.AddR,'Enable','on')
    set(handles.text10,'Enable','on')
    set(handles.edit4,'Enable','on')
    set(handles.text4,'Enable','on')
    set(handles.MaskThre,'Enable','on')
    set(handles.AddROIs,'String','Add all ROIs')
end



function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
