function varargout = AbortedHeader(varargin)
% Brunno Machado de Campos
% University of Campinas, 2016
%
% Copyright (c) 2016, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @AbortedHeader_OpeningFcn, ...
                   'gui_OutputFcn',  @AbortedHeader_OutputFcn, ...
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

function AbortedHeader_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = AbortedHeader_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function AddB_Callback(hObject, eventdata, handles)
global struI ImgO nVolReal

set(handles.text12,'String','Thinking...')
drawnow

[fileFunc,pathFunc] = uigetfile({'*.nii','NIfTI files'},...
    'Select the aborted functional image','MultiSelect','off');

set(handles.text12,'String','...')
set(handles.text12,'ForegroundColor',[1 0 0])
drawnow

ImgO = fullfile(pathFunc,fileFunc);
struI = nifti(ImgO);

matrixF = struI.dat.dim(1:3);
numofVoluF = struI.dat.dim(4);
VoxSizeF = struI.hdr.pixdim(2:4);
VoxS_Res = [num2str(VoxSizeF(1,1)),'x',num2str(VoxSizeF(1,2)),'x',num2str(VoxSizeF(1,3))];
MatS_Res = [num2str(matrixF(1,1)),'x',num2str(matrixF(1,2)),'x',num2str(matrixF(1,3))];
TRv = struI.timing.tspace;

fO=dir(ImgO);
BitO = fO.bytes;
struItmp = struI;
struItmp.dat.dim = [struI.dat.dim(1) struI.dat.dim(2) struI.dat.dim(3) 1];
Vol1 = struItmp.dat(:,:,:,1);
Stru2 = struI;

Stru2.dat.fname = [pathFunc filesep 'tempImg.nii'];
Stru2.dat.dim = Stru2.dat.dim(1:3);
Stru2.dat(:,:,:) = squeeze(Vol1);
create(Stru2);

f=dir([pathFunc filesep 'tempImg.nii']);
delete([pathFunc filesep 'tempImg.nii']);

Bit1 = f.bytes;
BitTot = Bit1.*numofVoluF;
BitPropo = BitO/BitTot;
nVolReal = floor(BitPropo.*numofVoluF);

duration = TRv.*numofVoluF;
duration = duration/60;
Hdur = floor(duration/60);
Mdur = floor(duration-(Hdur*60));
Sdur = (duration-floor(duration))*60;
TimeDur = [num2str(Hdur) ':' num2str(Mdur) ':' num2str(Sdur)];

set(handles.edit1,'String',MatS_Res);
set(handles.edit2,'String',numofVoluF);
set(handles.edit3,'String',VoxS_Res);
set(handles.edit5,'String',TRv);
set(handles.edit4,'String',TimeDur);

durationR = TRv.*nVolReal;
durationR = durationR/60;
HdurR = floor(durationR/60);
MdurR = floor(durationR-(HdurR*60));
SdurR = (durationR-floor(durationR))*60;
TimeDurR = [num2str(HdurR) ':' num2str(MdurR) ':' num2str(SdurR)];

set(handles.edit6,'String',MatS_Res);
set(handles.edit7,'String',nVolReal);
set(handles.edit8,'String',VoxS_Res);
set(handles.edit10,'String',TRv);
set(handles.edit9,'String',TimeDurR);

set(handles.text12,'ForegroundColor',[0 0 0])
set(handles.text12,'String',ImgO)
drawnow

function RunB_Callback(hObject, eventdata, handles)
global struI ImgO nVolReal

set(handles.text13,'String','Processing...')
drawnow

Stru3 = struI;
try
    Stru3.dat.dim = [Stru3.dat.dim(1:3) nVolReal];
    matF = Stru3.dat(:,:,:,:);
catch
    Stru3.dat.dim = [Stru3.dat.dim(1:3) nVolReal+1];
    matF = Stru3.dat(:,:,:,:);
end

Stru3.dat.fname = [ImgO(1:end-4) '_MOD_' get(handles.edit7,'String') 'vols.nii'];

Stru3.dat(:,:,:,:) = matF(:,:,:,1:str2num(get(handles.edit7,'String')));
create(Stru3);

set(handles.text13,'String','Done!')
drawnow


function edit10_Callback(hObject, eventdata, handles)

function edit10_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
