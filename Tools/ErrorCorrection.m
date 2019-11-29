function varargout = ErrorCorrection(varargin)

%Brunno Machado de Campos

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ErrorCorrection_OpeningFcn, ...
                   'gui_OutputFcn',  @ErrorCorrection_OutputFcn, ...
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

function ErrorCorrection_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = ErrorCorrection_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function erroex_Callback(hObject, eventdata, handles)
error_example

function Add_Callback(hObject, eventdata, handles)
global fname path
 
[fname,path] = uigetfile('*.nii','Add the image','MultiSelect','on');

if ~iscell(fname)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    fname = {fname};
end
fname = sort(fname);

for yt = 1:size(fname,2)
    tmpDir = dir([path fname{yt}]);
    bitS(yt) = tmpDir.bytes;
end
sxs = unique(bitS);
nDiftypes = size(sxs,2);
if ~isequal(nDiftypes,1)
    warndlg({sprintf('There are %d distinct protocol types among the images that you added.',nDiftypes);...
      sprintf('Check the size of your files, it is the easiest way to identify the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
end
clear bitS

verif = nifti([path,fname{1}]);
set(handles.addtext,'String',sprintf('%d image(s) added!',size(fname,2)))

if isequal(size(verif.dat.dim,2),3)
    numDy = '3D image';
    set(handles.checkbox3,'Value',0)
    set(handles.checkbox3,'Enable','off')
    set(handles.checkbox4,'Value',0)
    set(handles.checkbox4,'Enable','off')
    set(handles.checkbox5,'Value',0)
    set(handles.checkbox5,'Enable','off')
else
    numDy = verif.dat.dim(4);
    set(handles.checkbox3,'Enable','on')
    set(handles.checkbox4,'Enable','on')
    set(handles.checkbox5,'Enable','on')
    set(handles.toVolpre,'String',num2str(numDy))
    set(handles.toVolrem,'String',num2str(numDy))
end

set(handles.text20,'String',num2str(numDy))
set(handles.text22,'String',num2str(verif.dat.dim))
set(handles.text40,'String',num2str(verif.dat.dim(3)))
set(handles.text21,'String',num2str(verif.hdr.pixdim(2:4)))


function nslices_Callback(hObject, eventdata, handles)

function nslices_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function runb_Callback(hObject, eventdata, handles)
global fname path
disp('-------------------------------------------------------------------')

if isempty(get(handles.addtext,'String'))
    warndlg('First, you need to add a image file!','Ops!')
    return
end

for i = 1:size(fname,2)
    fprintf('Image %d of %d \n',i, size(fname,2));
    
    Naa  = nifti([path,fname{i}]);

    if isequal(get(handles.checkbox1,'Value'),1)
        nslice = str2double(get(handles.nslices,'String'));
        if isequal(nslice,0)
            warndlg('The value settled to number of mirrored slices is Zero','Attention!')
        else
            fprintf('Correcting slices...\n')
            if isequal(size(Naa.dat.dim,2),3) % 3D file
               xtrct = Naa.dat(:,:,1:nslice);
               nslc_tot = Naa.dat.dim(3); % total number of slices
               Y3D = Naa.dat(:,:,:);
               Naa.descrip = sprintf('%d slices put on top',nslice);
               Naa.dat.fname = [Naa.dat.fname(1:end-4) sprintf('_%dslc.nii',nslice)];
               Naa.dat(:,:,1:(nslc_tot-nslice)) = Y3D(:,:,(1+nslice):nslc_tot);
               Naa.dat(:,:,(nslc_tot-nslice+1):nslc_tot) = xtrct;

               create(Naa)

            else if isequal(size(Naa.dat.dim,2),4) % 4D file
                 xtrct = Naa.dat(:,:,1:nslice,:);
                 nslc_tot = Naa.dat.dim(3); % total number of slices
                 Y4D = Naa.dat(:,:,:,:);
                 Naa.descrip = sprintf('%d slices put on top',nslice);
                 Naa.dat.fname = [Naa.dat.fname(1:end-4) sprintf('_%dslc.nii',nslice)];
                 Naa.dat(:,:,1:(nslc_tot-nslice),:) = Y4D(:,:,(1+nslice):nslc_tot,:);
                 Naa.dat(:,:,(nslc_tot-nslice+1):nslc_tot,:) = xtrct;
                 create(Naa)
                end
            end
            disp('Done!')
        end
    end


    if isequal(get(handles.checkbox4,'Value'),1)
       if isequal(str2double(get(handles.SpecVol,'Value')),0)
          warndlg('The value settled to Volume number is Zero','Attention!')
       else
          fprintf('Removing volume...\n')
          espvol = str2double(get(handles.SpecVol, 'String'));
          newSize = Naa.dat.dim(4)-1;
          Y4D   = Naa.dat(:,:,:,:);
          Naa.dat.dim   = [Naa.dat.dim(1:3) newSize];
          Naa.descrip = sprintf('The original %dº volume was removed',espvol);
          Naa.dat.fname = [Naa.dat.fname(1:end-4) sprintf('_%d_rem.nii',espvol)];
          Y4D1 = Y4D(:,:,:,1:espvol-1);
          Y4D2 = Y4D(:,:,:,espvol+1:newSize+1);
          Naa.dat(:,:,:,1:espvol-1) = Y4D1;
          Naa.dat(:,:,:,espvol:newSize) = Y4D2;
          create(Naa)
          end
       disp('Done!')
    end

    if isequal(get(handles.checkbox3,'Value'),1)
       if isequal(str2double(get(handles.fromVolpre,'Value')),0) || isequal(str2double(get(handles.toVolpre,'Value')),0) ...
               || str2double(get(handles.fromVolpre,'Value'))>str2double(get(handles.toVolpre,'Value'))
          warndlg('The value settled to Volume interval is invalid','Attention!')
       else
          fprintf('Editing 4D image...\n')

          nPres = str2double(get(handles.toVolpre,'String'))-str2double(get(handles.fromVolpre,'String'));
          Y4D   = Naa.dat(:,:,:,str2double(get(handles.fromVolpre,'String')):str2double(get(handles.toVolpre,'String')));
          Naa.dat.dim   = [Naa.dat.dim(1:3) nPres+1];
          Naa.descrip = sprintf('Just %d volumes from the original images was preserved',nPres+1);
          Naa.dat.fname = [Naa.dat.fname(1:end-4) sprintf('_%d_pres.nii',nPres)];
          Naa.dat(:,:,:,1:nPres+1) = Y4D;
          create(Naa)
       end
       disp('Done!')
    end

    if isequal(get(handles.checkbox5,'Value'),1)
       if isequal(str2double(get(handles.fromVrem,'Value')),0) || isequal(str2double(get(handles.toVolrem,'Value')),0) ...
               || str2double(get(handles.fromVrem,'Value'))>str2double(get(handles.toVolrem,'Value'))
          warndlg('The value settled to Volume interval is invalid','Attention!')
       else
          fprintf('Editing 4D image...\n')

          nPres = str2double(get(handles.toVolrem,'String'))-str2double(get(handles.fromVrem,'String'))+1;
          SizeF = Naa.dat.dim(4)-nPres;
          Y4D1   = Naa.dat(:,:,:,1:str2double(get(handles.fromVrem,'String'))-1);
          Y4D2 = Naa.dat(:,:,:,str2double(get(handles.toVolrem,'String'))+1:end);

          Naa.dat.dim   = [Naa.dat.dim(1:3) Naa.dat.dim(4)-nPres];
          Naa.descrip = sprintf('%d volumes was removed from the original image.',nPres);
          Naa.dat.fname = [Naa.dat.fname(1:end-4) sprintf('_%d_remov.nii',nPres)];
          Naa.dat(:,:,:,1:str2double(get(handles.fromVrem,'String'))-1) = Y4D1;
          Naa.dat(:,:,:,str2double(get(handles.fromVrem,'String')):SizeF) = Y4D2;

          create(Naa)
       end
       disp('Done!')
    end
    clear Naa
end
disp('All Done!')

set(handles.addtext,'String','')

disp('-------------------------------------------------------------------')

function checkbox1_Callback(hObject, eventdata, handles)
    set(handles.checkbox3, 'Value' ,0)
    set(handles.checkbox4, 'Value' ,0)
    set(handles.checkbox5, 'Value' ,0)
    
    if isequal(get(handles.checkbox1,'Value'),1)
        set(handles.text28,'Enable','on')
        set(handles.erroex,'Enable','on')
        set(handles.nslices,'Enable','on')
        set(handles.runb,'Enable','on')
        
        set(handles.text32,'Enable','off')
        set(handles.SpecVol,'Enable','off')
        set(handles.text33,'Enable','off')
        set(handles.text34,'Enable','off')
        set(handles.fromVrem,'Enable','off')
        set(handles.toVolrem,'Enable','off')
        set(handles.text35,'Enable','off')
        set(handles.text36,'Enable','off')
        set(handles.fromVolpre,'Enable','off')
        set(handles.toVolpre,'Enable','off')
        set(handles.text37,'Enable','off')
        set(handles.text38,'Enable','off')
        set(handles.text41,'Enable','off')
        set(handles.text42,'Enable','off')
        set(handles.text43,'Enable','off')
        set(handles.text44,'Enable','off')
        
        
    else
        set(handles.text28,'Enable','off')
        set(handles.erroex,'Enable','off')
        set(handles.nslices,'Enable','off')
        set(handles.runb,'Enable','off')

    end

function checkbox3_Callback(hObject, eventdata, handles)
    set(handles.checkbox1, 'Value' ,0)
    set(handles.checkbox4, 'Value' ,0)
    set(handles.checkbox5, 'Value' ,0)
    
    if isequal(get(handles.checkbox3,'Value'),1)
        set(handles.text36,'Enable','on')
        set(handles.fromVolpre,'Enable','on')
        set(handles.toVolpre,'Enable','on')
        set(handles.text37,'Enable','on')
        set(handles.text38,'Enable','on')
        set(handles.runb,'Enable','on')
        set(handles.text41,'Enable','on')
        set(handles.text42,'Enable','on')
        
        set(handles.text33,'Enable','off')
        set(handles.text34,'Enable','off')
        set(handles.fromVrem,'Enable','off')
        set(handles.toVolrem,'Enable','off')
        set(handles.text35,'Enable','off')
        set(handles.text32,'Enable','off')
        set(handles.SpecVol,'Enable','off')
        set(handles.text28,'Enable','off')
        set(handles.erroex,'Enable','off')
        set(handles.nslices,'Enable','off')
        set(handles.text43,'Enable','off')
        set(handles.text44,'Enable','off')
        
    else
       set(handles.text36,'Enable','off')
       set(handles.fromVolpre,'Enable','off')
       set(handles.toVolpre,'Enable','off')
       set(handles.text37,'Enable','off')
       set(handles.text38,'Enable','off')
       set(handles.runb,'Enable','off')
       set(handles.text41,'Enable','off')
       set(handles.text42,'Enable','off')

    end
    


function checkbox4_Callback(hObject, eventdata, handles)
    set(handles.checkbox3, 'Value' ,0)
    set(handles.checkbox1, 'Value' ,0)
    set(handles.checkbox5, 'Value' ,0)

     if isequal(get(handles.checkbox4,'Value'),1)
        set(handles.text32,'Enable','on')
        set(handles.SpecVol,'Enable','on')
        set(handles.runb,'Enable','on')
            
        set(handles.text28,'Enable','of')
        set(handles.erroex,'Enable','of')
        set(handles.nslices,'Enable','of')
        set(handles.text33,'Enable','off')
        set(handles.text34,'Enable','off')
        set(handles.fromVrem,'Enable','off')
        set(handles.toVolrem,'Enable','off')
        set(handles.text35,'Enable','off')
        set(handles.text36,'Enable','off')
        set(handles.fromVolpre,'Enable','off')
        set(handles.toVolpre,'Enable','off')
        set(handles.text37,'Enable','off')
        set(handles.text38,'Enable','off')
        set(handles.text41,'Enable','off')
        set(handles.text42,'Enable','off')
        set(handles.text43,'Enable','off')
        set(handles.text44,'Enable','off')
    else
        set(handles.text32,'Enable','off')
        set(handles.SpecVol,'Enable','off')
        set(handles.runb,'Enable','off')
    end


function checkbox5_Callback(hObject, eventdata, handles)
    set(handles.checkbox3, 'Value' ,0)
    set(handles.checkbox4, 'Value' ,0)
    set(handles.checkbox1, 'Value' ,0)
    
    if isequal(get(handles.checkbox5,'Value'),1)
        set(handles.text33,'Enable','on')
        set(handles.text34,'Enable','on')
        set(handles.fromVrem,'Enable','on')
        set(handles.toVolrem,'Enable','on')
        set(handles.text35,'Enable','on')
        set(handles.runb,'Enable','on')
        set(handles.text43,'Enable','on')
        set(handles.text44,'Enable','on')
        
        set(handles.text41,'Enable','off')
        set(handles.text42,'Enable','off')
        set(handles.text36,'Enable','off')
        set(handles.text32,'Enable','off')
        set(handles.SpecVol,'Enable','off')
        set(handles.text28,'Enable','off')
        set(handles.erroex,'Enable','off')
        set(handles.nslices,'Enable','off')
        set(handles.fromVolpre,'Enable','off')
        set(handles.toVolpre,'Enable','off')
        set(handles.text37,'Enable','off')
        set(handles.text38,'Enable','off')
    else
        set(handles.text33,'Enable','off')
        set(handles.text34,'Enable','off')
        set(handles.fromVrem,'Enable','off')
        set(handles.toVolrem,'Enable','off')
        set(handles.text35,'Enable','off')
        set(handles.runb,'Enable','off')
        set(handles.text43,'Enable','off')
        set(handles.text44,'Enable','off')
    end

function fromVolpre_Callback(hObject, eventdata, handles)
Total3 = str2double(get(handles.toVolpre,'String'))-str2double(get(handles.fromVolpre,'String'))+1;
set(handles.text42,'String',num2str(Total3))
if isequal(str2double(get(handles.fromVolpre,'String')),0)
    set(handles.fromVolpre,'String','1')
end
Total3 = str2double(get(handles.toVolpre,'String'))-str2double(get(handles.fromVolpre,'String'))+1;
set(handles.text42,'String',num2str(Total3))

function fromVrem_Callback(hObject, eventdata, handles)
Total5 = str2double(get(handles.toVolrem,'String'))-str2double(get(handles.fromVrem,'String'))+1;
set(handles.text44,'String',num2str(Total5))
if isequal(str2double(get(handles.fromVolpre,'String')),0)
    set(handles.fromVolpre,'String','1')
end
Total5 = str2double(get(handles.toVolrem,'String'))-str2double(get(handles.fromVrem,'String'))+1;
set(handles.text42,'String',num2str(Total5))

function toVolrem_Callback(hObject, eventdata, handles)
Total5 = str2double(get(handles.toVolrem,'String'))-str2double(get(handles.fromVrem,'String'))+1;
set(handles.text44,'String',num2str(Total5))

function toVolpre_Callback(hObject, eventdata, handles)
Total3 = str2double(get(handles.toVolpre,'String'))-str2double(get(handles.fromVolpre,'String'))+1;
set(handles.text42,'String',num2str(Total3))




function fromVolpre_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toVolpre_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fromVrem_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toVolrem_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SpecVol_Callback(hObject, eventdata, handles)

function SpecVol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function VoluBeq_Callback(hObject, eventdata, handles)

function VoluBeq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nvolum_Callback(hObject, eventdata, handles)

function nvolum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
