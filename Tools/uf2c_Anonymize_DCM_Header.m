% Brunno Machado de Campos
% University of Campinas, 2018
%
% Copyright (c) 2018, Brunno Machado de Campos
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

function varargout = uf2c_Anonymize_DCM_Header(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_Anonymize_DCM_Header_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_Anonymize_DCM_Header_OutputFcn, ...
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


function uf2c_Anonymize_DCM_Header_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_Anonymize_DCM_Header_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddImg_Callback(hObject, eventdata, handles)
global fileFunc path2

[fileFunc,path2] = uigetfile({'*.*','DICOM files'},'Add all DICOM images','MultiSelect','on');

if ~iscell(fileFunc)
    fileFunc = {fileFunc};
end

set(handles.edit3,'String',[num2str(size(fileFunc,2)) ' images added!'])


function RunB_Callback(hObject, eventdata, handles)
global fileFunc path2

set(handles.text3,'String','Running....Wait!')

for i = 1:size(fileFunc,2)
    fprintf('Running...Wait!!\n')
    fprintf('Image %d of %d\n',i,size(fileFunc,2))
    
    metadata = alt_dicominfo([path2,fileFunc{i}]);

    FileName = '';

    if get(handles.checkbox9,'Value')
        FileName = [FileName '_' get(handles.edit1,'String')];
    end
    if get(handles.checkbox1,'Value')
        PatName = metadata.PatientName;
        SpacIdx = findstr(PatName,' ');
        Initi = [PatName(1) PatName(SpacIdx+1)];
        FileName = [FileName '_' Initi];
    end
    if get(handles.checkbox6,'Value')
        if isequal(metadata.PatientSex,'M')
            gender = 'Male';
        else
            gender = 'Female';
        end
        FileName = [FileName '_' gender];
    end
    if get(handles.checkbox7,'Value')
        dateF = [num2str(metadata.PatientBirth(end-1:end)),'-',...
            numb2month(str2num(metadata.PatientBirth(end-3:end-2)),...
            'short'),'-',num2str(metadata.PatientBirth(1:4))];

        dateFDT = datetime(dateF);

        dateStdy = [metadata.StudyDate(end-1:end),...
            numb2month(str2num(metadata.StudyDate(5:6)),'short'),...
            metadata.StudyDate(1:4)];
        dateStdyDT = datetime([dateStdy(1:2),'-',dateStdy(3:5),'-',dateStdy(6:9)]);
        nYears = num2str(round(years(dateStdyDT-dateFDT)));
        FileName = [FileName '_' nYears 'years'];
    end
    if get(handles.checkbox11,'Value')
        FileName = [FileName '_' metadata.PatientID];
    end
    if get(handles.checkbox3,'Value')
        Seridescr = metadata.SerieDescrip;
        FileName = [FileName '_' Seridescr];
    end
%     if get(handles.checkbox4,'Value')
%         FileName = [FileName '_' metadata.AcquisitionContrast];
%     end
%     if get(handles.checkbox4,'Value')
%         FileName = [FileName '_' [metadata.InstanceCreationDate(1:4) '-' metadata.InstanceCreationDate(5:6) '-' metadata.InstanceCreationDate(7:end)]];
%     end
%     if get(handles.checkbox5,'Value')
%         FileName = [FileName '_' metadata.InstitutionName];
%     end
    if get(handles.checkbox10,'Value')
        FileName = [FileName '_' get(handles.edit2,'String')];
    end
    if get(handles.checkbox12,'Value')
        FileName = [FileName '.dcm'];
    end
    if isequal(FileName(1),'_')
        FileName = FileName(2:end);
    end
    FileName = regexprep(FileName,' ','_');
    
    copyfile([path2,fileFunc{i}],[path2,FileName],'f')
    
    if get(handles.checkbox8,'Value')
        [~,fName,Ext] = fileparts(FileName);
        if get(handles.KPID,'Value')
            dicomanon([path2,FileName],[path2 fName '_Anon' Ext],...
                'keep',{'PatientBirthDate','InstitutionName','PatientID','PatientAge', 'PatientSex', 'StudyDescription'})
        else
            dicomanon([path2,FileName],[path2 fName '_Anon' Ext],...
                'keep',{'PatientBirthDate','InstitutionName','PatientAge', 'PatientSex', 'StudyDescription'})
        end            
            delete([path2,FileName])
    end

%     dicomwrite(X, [path2 FileName], metadata,'CreateMode','Copy','WritePrivate',true);
    clear X metadata pathA FileName AcquDate elapT Bdate PatName SpacIdx Initi
end
set(handles.text3,'String','Done!')
fprintf('Done! \n')


function checkbox8_Callback(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function checkbox3_Callback(hObject, eventdata, handles)

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)

function checkbox7_Callback(hObject, eventdata, handles)

function checkbox9_Callback(hObject, eventdata, handles)
if get(handles.checkbox9,'Value')
    set(handles.edit1,'Enable','on')
else
    set(handles.edit1,'Enable','off')
end

function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox10_Callback(hObject, eventdata, handles)
if get(handles.checkbox10,'Value')
    set(handles.edit2,'Enable','on')
else
    set(handles.edit2,'Enable','off')
end

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox11_Callback(hObject, eventdata, handles)

function checkbox12_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)


function KPID_Callback(hObject, eventdata, handles)
