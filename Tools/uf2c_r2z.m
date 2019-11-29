function varargout = uf2c_r2z(varargin)
% UF²C M-file for uf2c_r2z.fig
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
                   'gui_OpeningFcn', @uf2c_r2z_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_r2z_OutputFcn, ...
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

function uf2c_r2z_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_r2z_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function addb_Callback(hObject, eventdata, handles)

if get(handles.CorrSE,'Value')
    if isempty(get(handles.edit1,'String'))
        set(handles.outDLG,'String','You should to specify a Source Sample Size value')
        set(handles.text6,'String','Ops!')
        drawnow
        return
    end
end

[file,path] = uigetfile('*.*','Add NIfTI r-score maps or Matlab files with "double" values','MultiSelect','on');
set(handles.outDLG,'String','')
set(handles.text5,'String','')
set(handles.text6,'String','Running....Wait!')
drawnow

if ~iscell(file)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    file = {file};
end

for i = 1:size(file,2)
    file2{i} = [path file{i}];
end

[pa,name,ext] = fileparts(file2{1});

z = @(r) (0.5.*(log(1+r) - log(1-r)));
corrz = @(z) z./(1/sqrt(str2num(get(handles.edit1,'String'))-3));

for i = 1:size(file2,2)
    switch ext
        case '.nii'
            stru = nifti(file2{i});
            if isequal(size(stru.dat.dim),3)
                mat = stru.dat(:,:,:);
                mat(mat>0.99999999999999994) = 1;
                zmap(:,:,:) = z(mat(:,:,:));
                zmap(zmap==Inf) = 18.7150;

                if get(handles.CorrSE,'Value')
                   zmap = corrz(zmap);
                end
                
                stru.dat.fname = [path 'Z_Transf_' file{i}];
                stru.descrip = 'UF2C map. Z-score';
                stru.dat(:,:,:) = zmap;
                create (stru)
            else
                mat = stru.dat(:,:,:,:);
                mat(mat>0.99999999999999994) = 1;
                zmap(:,:,:,:) = z(mat(:,:,:,:));
                zmap(zmap==Inf) = 18.7150;
                
                if get(handles.CorrSE,'Value')
                     zmap = corrz(zmap);
                end

                stru.dat.fname = [path 'Z_Transf_' file{i}];
                stru.descrip = 'UF2C map. Z-score';
                stru.dat(:,:,:,:) = zmap;
                create (stru)
            end
        case '.mat'
            a = load(file2{i});
            value = fields(a);
            value = value{1};
            eval(sprintf('%s = a.%s;',value,value));
            Trans = eval(value);
            Trans(Trans>0.99999999999999994) = 1;
            Trans = z(Trans);
            Trans(Trans==Inf) = 18.7150;
            
            if get(handles.CorrSE,'Value')
                Trans = corrz(Trans);
            end
            eval(sprintf('%s = Trans;',value))
            save([path 'Z_Transf_' file{i}],sprintf('%s',value))
    end
end

set(handles.outDLG,'String','The transformed files were saved on the original files folders. All transformed files has the Z_Transf prefix')
set(handles.text5,'String','The variable name inside a .mat files (if added) was not modified.')
set(handles.text6,'String','Done!')
drawnow


function CorrSE_Callback(hObject, eventdata, handles)
if get(handles.CorrSE,'Value')
    set(handles.text7,'Enable','on')
    set(handles.edit1,'Enable','on')
else
    set(handles.text7,'Enable','off')
    set(handles.edit1,'Enable','off')
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


% --- Executes on button press in calcR.
function calcR_Callback(hObject, eventdata, handles)

set(handles.edit5,'String','');
set(handles.edit6,'String','');
set(handles.edit4,'String','');
set(handles.edit7,'String','');
drawnow

if isempty(get(handles.edit2,'String'))
    return
end
if isequal(get(handles.edit3,'String'),'0')
    set(handles.edit3,'String','')
end
% tval = @(r) r*sqrt((DF-2)/(1-r^2));
tval = @(r) r*sqrt((str2num(get(handles.edit3,'String'))-2)/(1-r^2));
z = @(r) (0.5.*(log(1+r) - log(1-r)));
corrz = @(z) z./(1/sqrt(str2num(get(handles.edit3,'String'))-3));

tvalor = tval(str2num(get(handles.edit2,'String')));

zmap = z(str2num(get(handles.edit2,'String')));
set(handles.edit4,'String',num2str(zmap));

if ~isempty(get(handles.edit3,'String'))
   SEv = 1/sqrt(str2num(get(handles.edit3,'String'))-3);
   zmapc = corrz(zmap);
   pvalor = 2*(1-tcdf(tvalor,str2num(get(handles.edit3,'String'))-2));
   set(handles.edit5,'String',num2str(SEv));
   set(handles.edit6,'String',num2str(zmapc));
   set(handles.edit7,'String',num2str(pvalor));
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



function edit6_Callback(hObject, eventdata, handles)
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
