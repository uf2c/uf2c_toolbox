function varargout = Maps_Intersection(varargin)

% Last Modified by GUIDE v2.5 31-May-2017 12:54:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Maps_Intersection_OpeningFcn, ...
                   'gui_OutputFcn',  @Maps_Intersection_OutputFcn, ...
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
% End initialization code

function Maps_Intersection_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = Maps_Intersection_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function VBM_Callback(hObject, eventdata, handles)
global file1o path1r filesFold

if get(handles.pair,'Value')
    [file1o,path1r] = uigetfile({'*.nii;*.img','IMAGES'},'Add the map with higher resolution');
    if ~isempty(file1o)
        set(handles.textfiles1,'String',file1o)
    end
else
   filesFold = uipickfiles('Output','cell','Prompt','Add all subjects directories','REFilter','\');
    if ~isempty(filesFold)
        set(handles.textfiles1,'String',[num2str(numel(filesFold)) ' folder(s) added'])
        drawnow
    end
end

function BOLD_Callback(hObject, eventdata, handles)
global file2o path2r

[file2o,path2r] = uigetfile({'*.nii;*.img','IMAGES'},'Add the map with lower or equal resolution');
if ~isempty(file2o)
    set(handles.textfiles2,'String',file2o)
end

function run_Callback(hObject, eventdata, handles)
global file1o path1r file2o path2r filesFold

spm_jobman('initcfg')

if get(handles.pair,'Value')
    loop = 1;
    pathRes = path1r;
    fide = fopen([pathRes 'Intersection results.txt'],'w+');
else
    loop = numel(filesFold);
    pathRes = [filesFold{1} filesep];
    fide = fopen([pathRes 'Intersection results.txt'],'w+');
end
set(handles.text2,'String','')
set(handles.text12,'String', 'Running.....Wait!');
drawnow

for i = 1:loop
    
    clear matlabbatch
    Invert = 0;
    if get(handles.pair,'Value')
        file1r = fullfile(path1r,file1o);
        file2r = fullfile(path2r,file2o);
    else
        a = dir(filesFold{i});
        idx = 1;
        for j = 3:numel(a)
            [p,na,ext] = fileparts(a(j,1).name);
            if isequal(ext,'.nii') || isequal(ext,'.img')
                files{idx,1} = a(j,1).name;
                idx = 2;
            end
        end
        file1r = [filesFold{i} filesep files{1}];
        path1r = [filesFold{i} filesep];
        file2r = [filesFold{i} filesep files{2}];
        path2r = [filesFold{i} filesep];
    end

    Stru1 = nifti(file1r);
    tmp1pix = Stru1.hdr.pixdim(2:4);
    tmp1dim = Stru1.dat.dim(1:3);

    Stru2 = nifti(file2r);
    tmp2pix = Stru2.hdr.pixdim(2:4);
    tmp2dim = Stru2.dat.dim(1:3);

    if prod(tmp2pix) < prod(tmp1pix)
        Invert = 1;
        file2 = file1r;
        file1 = file2r;

        Stru1 = nifti(file1);
        tmp1pix = Stru1.hdr.pixdim(2:4);
        tmp1dim = Stru1.dat.dim(1:3);

        Stru2 = nifti(file2);
        tmp2pix = Stru2.hdr.pixdim(2:4);
        tmp2dim = Stru2.dat.dim(1:3);
    else
        file1 = file1r;
        file2 = file2r;
    end

    if ~isequal(tmp1pix,tmp2pix) || ~isequal(tmp1dim,tmp2dim)
        
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {file1};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {file2};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        spm_jobman('run',matlabbatch)
        clear matlabbatch

        if get(handles.pair,'Value')
            if Invert
                file2 = [path1r filesep 'r' file1o];
                Stru2 = nifti(file2);
            else
                file2 = [path2r filesep 'r' file2o];
                Stru2 = nifti(file2);
            end
        else
            if Invert
                file2 = [path1r filesep 'r' files{1}];
                Stru2 = nifti(file2);
            else
                file2 = [path2r filesep 'r' files{2}];
                Stru2 = nifti(file2);
            end
        end
    end

    mat1 = Stru1.dat(:,:,:);
    mat1(mat1>str2num(get(handles.edit1,'String'))) = 1;
    mat1(mat1<=str2num(get(handles.edit1,'String'))) = 0;
    mat1(isnan(mat1)) = 0;
    TotVox1 = sum(sum(sum(mat1)));

    mat2 = Stru2.dat(:,:,:);
    mat2(mat2>str2num(get(handles.edit2,'String'))) = 1;
    mat2(mat2<=str2num(get(handles.edit2,'String'))) = 0;
    mat2(isnan(mat2)) = 0;
    TotVox2 = sum(sum(sum(mat2)));

    IntercMap = mat1.*mat2;
    TotVoxInter = sum(sum(sum(IntercMap)));

    Perc2file1 = 100*(TotVoxInter/TotVox1);
    Perc2file2 = 100*(TotVoxInter/TotVox2);

    if get(handles.multf,'Value')
        fprintf(fide,'Folder:\t %s  \r\n',filesFold{i});
    end
    
    [x1,file1,z1] = fileparts(file1);
    [x2,file2,z2] = fileparts(file2);
    
    if isnan(Perc2file1)
        Perc2file1 = 0;
    end
    if isnan(Perc2file2)
        Perc2file2 = 0;
    end

    fprintf(fide,'\t Voxels of %s which match %s: \t %.2f%%\r\n',file1,file2,Perc2file1);
    fprintf(fide,'\t Voxels of %s which match %s: \t %.2f%%\r\n',file2,file1,Perc2file2);
    fprintf(fide,'\t Number of voxel in %s: \t %d\r\n',file1,TotVox1);
    fprintf(fide,'\t Number of voxel in %s: \t %d\r\n',file2,TotVox2);
    fprintf(fide,'\t Total intersectioned voxels: \t %d\r\n\r\n',TotVoxInter);
    
%     delete([x2 filesep file2 z2])
end
fclose('all');
set(handles.text2,'String',['Results saved on: ' pathRes 'Intersection_results.txt'])
set(handles.text12,'String', 'Done!!');

function pair_Callback(hObject, eventdata, handles)
if get(handles.pair,'Value')
    set(handles.multf,'Value',0)
    set(handles.BOLD,'Enable','on')
    set(handles.VBM,'String','First Image')
else
    set(handles.multf,'Value',1)
    set(handles.BOLD,'Enable','off')
    set(handles.VBM,'String','Add Folders')
end
function multf_Callback(hObject, eventdata, handles)
if get(handles.multf,'Value')
    set(handles.pair,'Value',0)
    set(handles.BOLD,'Enable','off')
    set(handles.VBM,'String','Add Folders')
    set(handles.VBM,'TooltipString','Add Folders with only the pair of images. One folder per subject/case')
else
    set(handles.pair,'Value',1)
    set(handles.BOLD,'Enable','on')
    set(handles.VBM,'String','First Image')
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
