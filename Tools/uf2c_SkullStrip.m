function varargout = uf2c_SkullStrip(varargin)
% UF2C_SKULLSTRIP MATLAB code for uf2c_SkullStrip.fig
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_SkullStrip_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_SkullStrip_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before uf2c_SkullStrip is made visible.
function uf2c_SkullStrip_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);
function varargout = uf2c_SkullStrip_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Yesseg_Callback(hObject, eventdata, handles)
set(handles.Noseg, 'Value', 0)

if get(handles.Yesseg, 'Value') == 0
    set(handles.Noseg, 'Value', 1)
end

if get(handles.Yesseg,'Value') ==1
    set(handles.Addimages, 'Visible', 'off')
    set(handles.Addseg, 'Visible', 'on')
    set(handles.RunYseg, 'Visible','on')
    set(handles.run, 'Visible','off')
    set(handles.addrawimag, 'Enable','on')
    set(handles.text5, 'Enable','off')
    set(handles.x, 'Enable','off')
    set(handles.y, 'Enable','off')
    set(handles.z, 'Enable','off')
    set(handles.thres, 'String','0.9')

end
if get(handles.Yesseg,'Value') ==0
    set(handles.Addimages, 'Visible', 'on')
    set(handles.Addseg, 'Visible', 'off')
    set(handles.RunYseg, 'Visible','off')
    set(handles.run, 'Visible','on')
    set(handles.addrawimag, 'Enable','off')
    set(handles.text5, 'Enable','on')
    set(handles.x, 'Enable','on')
    set(handles.y, 'Enable','on')
    set(handles.z, 'Enable','on')
    set(handles.thres, 'String','0.32')

end

% --- Executes on button press in Noseg.
function Noseg_Callback(hObject, eventdata, handles)
set(handles.Yesseg, 'Value', 0)
if get(handles.Noseg, 'Value') == 0
    set(handles.Yesseg, 'Value', 1)
end
if get(handles.Yesseg,'Value') ==1
    set(handles.Addimages, 'Visible', 'off')
    set(handles.Addseg, 'Visible', 'on')
    set(handles.RunYseg, 'Visible','on')
    set(handles.run, 'Visible','off')
    set(handles.addrawimag, 'Enable','on')
    set(handles.text5, 'Enable','off')
    set(handles.x, 'Enable','off')
    set(handles.y, 'Enable','off')
    set(handles.z, 'Enable','off')
    set(handles.thres, 'String','0.9')

end
if get(handles.Yesseg,'Value') ==0
    set(handles.Addimages, 'Visible', 'on')
    set(handles.Addseg, 'Visible', 'off')
    set(handles.RunYseg, 'Visible','off')
    set(handles.run, 'Visible','on')
    set(handles.addrawimag, 'Enable','off')
    set(handles.text5, 'Enable','on')
    set(handles.x, 'Enable','on')
    set(handles.y, 'Enable','on')
    set(handles.z, 'Enable','on')
    set(handles.thres, 'String','0.32')

end

function thres_Callback(hObject, eventdata, handles)

function thres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Addseg_Callback(hObject, eventdata, handles)
global file1 path1

[file1, path1] = uigetfile('*.nii','Select the segmented images','MultiSelect','on');
nf = cellstr(file1);
nf = size(nf);
nf = nf(1,2);
set(handles.showfile, 'String', sprintf('%d files addes',nf));

function Addimages_Callback(hObject, eventdata, handles)
global fileS pathS

[fileS, pathS] = uigetfile('*.*','Select the raw image');
set(handles.showfile, 'String', fileS);

function run_Callback(hObject, eventdata, handles)
global fileS  pathS 

set(handles.RunningTxt, 'String','Wait...Running....');

a = str2double(get(handles.x, 'String'));
b = str2double(get(handles.y,'String'));
c = str2double(get(handles.z,'String'));

aux2 = [a b c];
            
            clear matlabbatch 
            matlabbatch{1}.spm.tools.preproc8.channel.vols = {[pathS fileS]};
            matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
            matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
            matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
            
            versi = spm('Ver','',1);
    
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,1']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,1']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
            
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,2']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,2']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
            
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,3']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,3']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
            
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,4']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,4']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
            
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,5']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,5']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];             
            
            if isequal(versi,'SPM8')
                matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('dir') filesep 'toolbox' filesep 'Seg' filesep 'TPM.nii,6']};           
            elseif isequal(versi,'SPM12b')
                matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('dir') filesep 'tpm' filesep 'TPM.nii,6']};
            end
            
            matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [1 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
            matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
            matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
            matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];
            save('batch.mat','matlabbatch')
            spm_jobman('run',matlabbatch)
            
            file2 = {[pathS 'c1' fileS];[pathS 'c2' fileS];[pathS 'c3' fileS];[pathS 'c4' fileS];[pathS 'c5' fileS];[pathS 'c6' fileS]};
              
            aux = str2double(get(handles.thres,'String'));
            
            clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input          = file2;
            matlabbatch{1}.spm.util.imcalc.output         = 'mask.nii';
            matlabbatch{1}.spm.util.imcalc.outdir         = {''};
            matlabbatch{1}.spm.util.imcalc.expression     = sprintf('((i1>0.5)+(i2>0.5)+(i3>0.5))-((i4)+(i5)+(i6))');
            matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
            spm_jobman('run',matlabbatch)
            
            file2 = {[pathS 'c4' fileS];[pathS 'c5' fileS];[pathS 'c6' fileS]};
            clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input          = file2;
            matlabbatch{1}.spm.util.imcalc.output         = 'mask.nii';
            matlabbatch{1}.spm.util.imcalc.outdir         = {''};
            matlabbatch{1}.spm.util.imcalc.expression     = 'i1+i2+i3';
            matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
            spm_jobman('run',matlabbatch)
            
            
            
            clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input          = {[pathS fileS];[pathS 'mask.nii']};
            matlabbatch{1}.spm.util.imcalc.output         = 'Skull Striped.nii';
            matlabbatch{1}.spm.util.imcalc.outdir         = {''};
            matlabbatch{1}.spm.util.imcalc.expression     = 'i1.*i2';
            matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;

            spm_jobman('run',matlabbatch)
            set(handles.RunningTxt, 'String','Done!');



function x_Callback(hObject, eventdata, handles)

function x_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_Callback(hObject, eventdata, handles)

function y_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function z_Callback(hObject, eventdata, handles)

function z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RunYseg_Callback(hObject, eventdata, handles)

global file1  path1 file2 path2

aux = str2double(get(handles.thres,'String'));
            
            cd (path1)

            
file1 = cellstr(file1);
            file1 = transpose(file1);
            num = numel([file1]);    % gera variavel com o tamanho da matriz
            
            A{1,1} = 'i1';
            
            for j=2:num             % cria matriz A do tamanho de num
                A{1,j} =  sprintf('+i%d',j);
            end
            
            EQUA = cat(2,A{:});     %Concatena valores das celulas da matriz A(1,num)
           
            set(handles.RunningTxt, 'String','Wait...Running....');


            matlabbatch{1}.spm.util.imcalc.input          = file1;
            matlabbatch{1}.spm.util.imcalc.output         = 'mask.nii';
            matlabbatch{1}.spm.util.imcalc.outdir         = {''};
            matlabbatch{1}.spm.util.imcalc.expression     = sprintf('%s>%.1f',EQUA,aux);
            matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;

            
            
            matlabbatch{2}.spm.util.imcalc.input          = {fullfile(path2,file2);fullfile(pwd,'mask.nii')};
            matlabbatch{2}.spm.util.imcalc.output         = 'Skull Striped.nii';
            matlabbatch{2}.spm.util.imcalc.outdir         = {''};
            matlabbatch{2}.spm.util.imcalc.expression     = 'i1.*i2';
            matlabbatch{2}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{2}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{2}.spm.util.imcalc.options.interp = 1;
            matlabbatch{2}.spm.util.imcalc.options.dtype  = 4;

            spm_jobman('run',matlabbatch)
            set(handles.RunningTxt, 'String','Done!');
            
            
            

% --- Executes on button press in addrawimag.
function addrawimag_Callback(hObject, eventdata, handles)
global file2 path2
[file2, path2] = uigetfile('*.*','Select the raw images');
set(handles.rawimag,'String',file2)
