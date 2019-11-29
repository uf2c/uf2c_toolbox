function varargout = uf2c_InterpGUI(varargin)
% Bruno Machado de Campos
% University of Campinas

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_InterpGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_InterpGUI_OutputFcn, ...
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

function uf2c_InterpGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_InterpGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddBtm_Callback(hObject, eventdata, handles)
global fileStr pathStr
[fileStr,pathStr] = uigetfile({'*.nii' ,'*.nii (NIfTI)'},...
    'Select all target files' ,'MultiSelect','on');

if ~iscell(fileStr)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    fileStr = {fileStr};
end
set(handles.AddTxt,'String','Image(s) added!')

function RunBtm_Callback(hObject, eventdata, handles)
global fileStr pathStr fileEx pathEx

set(handles.text4,'String','Running...')
drawnow
if isequal(get(handles.popupmenu1,'Value'),1)
    pathDD = which('uf2c');
    pathDD = [pathDD(1:end-6) 'Image_Examples' filesep 'MNI_StdRes_Mask.nii'];

    for i = 1:size(fileStr,2)
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {
                                        pathDD
                                        [pathStr filesep fileStr{i}]
                                        };
        matlabbatch{1}.spm.util.imcalc.output = [pathStr filesep fileStr{i}(1:end-4) '_Interp.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {pathStr};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        if get(handles.radiobutton1,'Value')
            matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        else
            matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        end
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)
    end
else
    
    pathDD = [pathEx filesep fileEx];
    Struex = nifti(pathDD);
    matex = Struex.dat(:,:,:);
    matex(:,:,:) = 1;
    tmpS = Struex;
    tmpS.dat.fname = [Struex.dat.fname(1:end-4) '_TempMask.nii'];
    tmpS.dat(:,:,:) = matex;
    create(tmpS)
    
    pathDD = [Struex.dat.fname(1:end-4) '_TempMask.nii'];
    
    for i = 1:size(fileStr,2)
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {
                                        pathDD
                                        [pathStr filesep fileStr{i}]
                                        };
        matlabbatch{1}.spm.util.imcalc.output = [pathStr filesep fileStr{i}(1:end-4) '_Interp.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {pathStr};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        if get(handles.radiobutton1,'Value')
            matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        else
            matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        end
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)
    end
    delete([Struex.dat.fname(1:end-4) '_TempMask.nii'])
end
set(handles.text4,'String','Done!')    
drawnow


function popupmenu1_Callback(hObject, eventdata, handles)
if isequal(eventdata.Source.Value,1)
    set(handles.CustB,'Enable','off')
    set(handles.text7,'Enable','off')
else
    set(handles.CustB,'Enable','on')
    set(handles.text7,'Enable','on')
end

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CustB_Callback(hObject, eventdata, handles)
global fileEx pathEx

[fileEx,pathEx] = uigetfile({'*.nii' ,'*.nii (NIfTI)'},...
    'Select all target files' ,'MultiSelect','off');

set(handles.text7,'String','Image added!')


function radiobutton1_Callback(hObject, eventdata, handles)

if get(handles.radiobutton1,'Value')
    set(handles.radiobutton2,'Value',0)
else
    set(handles.radiobutton2,'Value',1)
end
    
    
function radiobutton2_Callback(hObject, eventdata, handles)
if get(handles.radiobutton2,'Value')
    set(handles.radiobutton1,'Value',0)
else
    set(handles.radiobutton1,'Value',1)
end
