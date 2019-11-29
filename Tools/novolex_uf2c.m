function varargout = novolex_uf2c(varargin)
% UF²C M-file for novolex_uf2c.fig
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
                   'gui_OpeningFcn', @novolex_uf2c_OpeningFcn, ...
                   'gui_OutputFcn',  @novolex_uf2c_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

warning('off','all')

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function novolex_uf2c_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = novolex_uf2c_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddFunc_Callback(hObject, eventdata, handles)
global fileFunc pathFunc numofVoluF nOFsubjects matrixF VoxSizeF

set(handles.checkReg,'Enable','on')

try
    clear('regVets');
end
try
    delete([tmpDIR 'additional_Reg.mat']);
end

[fileFunc,pathFunc] = uigetfile({'*.nii','NIfTI files'},'Select all the functional images','MultiSelect','on');

if ~isequal(fileFunc,0)
    if ~iscell(fileFunc)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileFunc = {fileFunc};
    end
    fileFunc = sort(fileFunc); % SORT FILES IN THE ALPHABETIC ORDER
    nOFsubjects = size(fileFunc,2);
    for yt = 1:nOFsubjects
        tmpDir = dir([pathFunc fileFunc{yt}]);
        bitS(yt) = tmpDir.bytes;
    end
    sxs = unique(bitS);
    nDiftypes = size(sxs,2);
    if ~isequal(nDiftypes,1)
        warndlg({sprintf('There are %d distinct protocol types among the functional images that you added.',nDiftypes);...
          sprintf('Check the size of your files, it is the easiest way to identify the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
    end
    clear bitS

    set(handles.txtFunc,'String',sprintf('%d functional image(s) added',nOFsubjects))

    previewF = nifti([pathFunc fileFunc{1}]);
    matrixF = previewF.dat.dim(1:3);
    numofVoluF = previewF.dat.dim(4);
    VoxSizeF = previewF.hdr.pixdim(2:4);
    VoxS_Res = [num2str(VoxSizeF(1,1)),'x',num2str(VoxSizeF(1,2)),'x',num2str(VoxSizeF(1,3))];
    MatS_Res = [num2str(matrixF(1,1)),'x',num2str(matrixF(1,2)),'x',num2str(matrixF(1,3))];

    set(handles.numDyna,'String',numofVoluF)
    set(handles.SizeVoxFunc,'String',VoxS_Res)
    set(handles.SizeMatFunc,'String',MatS_Res)
    set(handles.LVR,'String',floor(numofVoluF/10))
end

function AddStru_Callback(hObject, eventdata, handles)
global fileStru pathStru

[fileStru,pathStru] = uigetfile({'*.nii','NIfTI files'},'Select all the Structural images','MultiSelect','on');

if ~isequal(fileStru,0)
    if ~iscell(fileStru)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileStru = {fileStru};
    end

    fileStru = sort(fileStru);
end

function Run_Callback(hObject, eventdata, handles)
global fileStru pathStru fileFunc pathFunc numofVoluF nrg tmpDIR

multiWaitbar('CloseAll'); %close the progressbar window

if ~isequal(size(fileStru,2),size(fileFunc,2))
    warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
      'You need to add one functional image for each strutural and vice versa.';...
      'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
      'Attention: process aborted!');
    return
end

set(handles.status,'String','Running....')
drawnow

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Process  started\n');
fprintf('==========================================\r\n');

multiWaitbar('Total Progress', 0, 'Color', 'g' );  %CHARGE WAIT BAR
multiWaitbar('Total Preprocessing Progress', 0, 'Color', 'y' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('1-Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

imgRR = getframe(Preproc_RMC);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Imaging Processing and Quantifications\r\n');

for yy = 1: size(fileFunc,2) % START OF SUBJECT LOOP
    
    multiWaitbar('Total Preprocessing Progress', (yy/size(fileFunc,2)).*0.95);
    
    fprintf('\r\n')
    fprintf('UF²C ===================== %s\n',spm('time'));
    if size(fileFunc{yy},2)<15
        fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy});
    else
        fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy}(1:14));
    end
    fprintf('================================================\r\n');
    
    try
        close(fig)
    end
    
    Funpre = nifti([pathFunc,fileFunc{yy}]);
    fvs = double(Funpre.hdr.pixdim(2:4)); %GET PIXEL SIZE (FUNCTIONAL)
    dirname = Funpre.dat.fname;
    [ax,dirname,cx] = fileparts(dirname);
    
    if numel(dirname)>30
        fprintf('Attention!\n');
        fprintf('The filename is long. A shorter version will be used.\r\n')
        dirname = dirname(1:30);
    end
    
    nfidx = 1;
    while isequal(exist([pathFunc,dirname],'dir'),7)
        dirname = [dirname '_' num2str(nfidx)];
        nfidx = nfidx + 1;
    end
    
    mkdir([pathFunc,dirname])
    fprintf(Totmovtxt,'%s \t',dirname);

    copyfile([pathStru,fileStru{yy}],[pathFunc,dirname,filesep,fileStru{yy}])
    
    file2 = fullfile(pathFunc,dirname,filesep,fileStru{yy});
    copyfile([pathFunc,fileFunc{yy}],[pathFunc,dirname,filesep,fileFunc{yy}])

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Removing Suprathreshold Volumes\r\n');
    fprintf('Step 1:');
    
    
    
    

    fprintf('Done!\r\n');

    fprintf('Step 3: ');

    idxs = 1;
    IncluSubj = [];
    
    fprintf('Subject %d Done! ====== %s\n\r',yy,spm('time'));
    
end

save([pathFunc,foldername,filesep,'StruTOT'],'StruTOT');

multiWaitbar('Total Filtering & Regression Process', 1);
multiWaitbar('Total Progress', 1);
multiWaitbar('CloseAll'); %close the progressbar window

set(handles.status,'String','Done!!!')

fclose('all');

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('All Done!\n');
fprintf('==========================================\r\n');


function checFunc_Callback(hObject, eventdata, handles)
global fileFunc

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Func.txt']);
end

subL = transpose(fileFunc);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Func.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileFunc{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Func.txt'])

function checStru_Callback(hObject, eventdata, handles)
global fileStru

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Stru.txt']);
end

subL = transpose(fileStru);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Stru.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileStru{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Stru.txt'])

function DVARthres_Callback(hObject, eventdata, handles)
function FDthres_Callback(hObject, eventdata, handles)

function FDthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DVARthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extGM_Callback(hObject, eventdata, handles)
if get(handles.extGM,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)

else
    if get(handles.extWM,'Value') || get(handles.extCSF,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)

    end
    if isequal(get(handles.extCSF,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extGM,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    end
end
    
function extWM_Callback(hObject, eventdata, handles)
if get(handles.extWM,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extCSF,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extCSF,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extWM,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    end

end

function extCSF_Callback(hObject, eventdata, handles)
if get(handles.extCSF,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extWM,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)

    end
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extCSF,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)

    end
end

function FDcheck_Callback(hObject, eventdata, handles)
if get(handles.FDcheck,'Value')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.extCSF,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.FDcheck,'Value',1)
        set(handles.text28,'Enable','on')
        set(handles.FDthres,'Enable','on')

    end
end

function LVR_Callback(hObject, eventdata, handles)

function LVR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
