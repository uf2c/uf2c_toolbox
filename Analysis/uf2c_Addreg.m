function varargout = uf2c_Addreg(varargin)
% UF²C M-file for FuncCon.fig
% UF²C - User frendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2013

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_Addreg_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_Addreg_OutputFcn, ...
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

function uf2c_Addreg_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_Addreg_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function addtregs_Callback(hObject, eventdata, handles)

function addtregs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SandCb_Callback(hObject, eventdata, handles)
global regVets tmpDIR

UF2Cdir = which('uf2c_FuncCon');
tmpDIR = [UF2Cdir(1:end-14) 'FC_tmp' filesep];
try
    save([tmpDIR 'additional_Reg'],'regVets')
catch
    warndlg('No additional regressors added!')
end
close(uf2c_Addreg)

function AddregB_Callback(hObject, eventdata, handles)
global numofVoluF nrg nOFsubjects regVets

[fileReg,pathReg] = uigetfile({'*.xlsx','XLSX files';'*.xls','XLS files'},'Select the XLSX (or XLS) file','MultiSelect','off');


[NUM,TXT,RAW] = xlsread([pathReg fileReg]);

regVets = zeros(size(RAW,1),size(RAW,2));

try
    for i = 1:size(regVets,2)
        for j = 1:size(regVets,1)
            if ischar(RAW{j,i})
                regVets(j,i) = str2num(RAW{j,i});
            else
                regVets(j,i) = RAW{j,i};
            end
        end
    end
catch
    warndlg('There is an error in your XLS file')
end

nOFsub = size(RAW,2);
if ~isequal(nOFsub,nOFsubjects)
    warndlg('The number of columns in you table do not match with the number of subjects previously added')
end

nOFReg = size(RAW,1)/numofVoluF;
nOFRegT = num2str(nOFReg);

if  ~isequal(round(nOFReg),nOFReg)
    warndlg('The number of lines (regressor series) is not multiple of the number of Volumes (functional series)')
    set(handles.nregTxt,'String','Error');
else
    set(handles.nregTxt,'String',nOFRegT);
    nrg = nOFReg;
end

verN = 1;
linN = 1;

for tr = 1:numofVoluF:size(regVets,1)
    for gh = 1:size(regVets,2)
        if  size(unique(regVets(tr:tr+numofVoluF-1,gh)),1)==1
            badReg(linN,:) = [verN gh];
            linN = linN+1;
        end
    end
    verN = verN+1;
end

if exist('badReg')
    for iu = 1: size(badReg,1)
        strErr(iu,1) = {sprintf('The regressor %d from subject %d seems to have all same values.',badReg(iu,1),badReg(iu,2))};
    end
    strErr(iu+2,1) = {'All listed regressors will be not considered!'};
warndlg(strErr);
end

function Clcregb_Callback(hObject, eventdata, handles)
global nrg regVets tmpDIR
UF2Cdir = which('FuncCon');
tmpDIR = [UF2Cdir(1:end-9) 'FC_tmp' filesep];

nrg = 0;

try
    clear('regVets');
end

try
    delete([tmpDIR 'additional_Reg.mat']);
end

set(handles.clctxt,'String','cleaned');
set(handles.nregTxt,'String','0');

function refND_Callback(hObject, eventdata, handles)
global nOFsubjects numofVoluF

set(handles.dynN,'String',numofVoluF);
set(handles.nsubTXT,'String',nOFsubjects);
set(handles.AddregB,'Enable','on')
    
function nsubTXT_Callback(hObject, eventdata, handles)

function clctxt_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nregTxt_Callback(hObject, eventdata, handles)

function nregTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function readInst_Callback(hObject, eventdata, handles)
msgbox({'1 - The regressors should to be added in a XLSX (or XLS) file.';...
    '2 - One columns for each subject included (see subject order list).';...
    '3 - The subject order and the regressor column should to match.';...
    '4 - All different regressors should to be included sequentialy in the same column.';...
    '5 - Your can add how many regressors do you want (but the same number for all subjecs)'},'Instructions');

function SubList_Callback(hObject, eventdata, handles)
global fileFunc tmpDIR

UF2Cdir = which('FuncCon');
tmpDIR = [UF2Cdir(1:end-9) 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List.txt']);
end

subL = transpose(fileFunc);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List.txt'],'w+'); % CHARGE OUTPUT LOG FILE

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileFunc{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List.txt'])

function nsubTXT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
