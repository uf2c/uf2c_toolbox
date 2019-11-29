function varargout = uf2c_Volumetry(varargin)
% UF²C M-file for uf2c.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2019
%
% Copyright (c) 2019, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @uf2c_Volumetry_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_Volumetry_OutputFcn, ...
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
function uf2c_Volumetry_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = uf2c_Volumetry_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function Addfiles_Callback(hObject, eventdata, handles)
global filevol pathvol

[filevol,pathvol] = uigetfile({'*.nii','*.nii (NIfTI)';'*.img','*.img (Analyze)'},...
    'Select all segmented images (grey matter, white matter and CSF) of all volumes', 'MultiSelect', 'on');

if ~iscell(filevol)
    filevol = {filevol};
    set(handles.textfiles, 'String', 'file added')
else
    set(handles.textfiles, 'String', 'files added')
end

AAA = size(filevol,2);
set(handles.nfiles, 'String', num2str(AAA))
nOfSub = AAA/3;

if get(handles.threefiles, 'Value')
     set(handles.textfiles, 'String', ['files added (' num2str(nOfSub) ' subjects)'])
end


function Run_Callback(hObject, eventdata, handles)
global filevol pathvol
   
    set(handles.text12,'String', 'Running.....Wait!');
    drawnow
    
    fide = fopen([pathvol filesep 'Volumetry_Results.txt'],'w+');
    AAA = size(filevol,2);
    nOfSub = AAA/3;
    
    for i = 1:nOfSub

        fprintf('Volume %d\n',i);
        
        C1map = [pathvol filesep filevol{i}];
        C2map = [pathvol filesep filevol{i+nOfSub}];
        C3map = [pathvol filesep filevol{i+(2*nOfSub)}];
        
        StruC1 = nifti(C1map);
        MatC1 = StruC1.dat(:,:,:);
        MatC1(MatC1<str2num(get(handles.thrs,'String'))) = 0;
        MatC1(MatC1>=str2num(get(handles.thrs,'String'))) = 1;
        
        VoxDim = StruC1.hdr.pixdim(2:4);
        
        StruC2 = nifti(C2map);
        MatC2 = StruC2.dat(:,:,:);
        MatC2(MatC2<str2num(get(handles.thrs,'String'))) = 0;
        MatC2(MatC2>=str2num(get(handles.thrs,'String'))) = 1;

        StruC3 = nifti(C3map);
        MatC3 = StruC3.dat(:,:,:);
        MatC3(MatC3<str2num(get(handles.thrs,'String'))) = 0;
        MatC3(MatC3>=str2num(get(handles.thrs,'String'))) = 1;

        
        TotMask = MatC1 + MatC2 + MatC3;
        TotVol = nnz(TotMask) * prod(VoxDim);
            TotVol_cm = TotVol/1000;
                TotVol_lts = TotVol/1000000;
                
        C1Vol = nnz(MatC1) * prod(VoxDim);
            C1Vol_cm = C1Vol/1000;
                C1Vol_lts = C1Vol/100000;
                
        C2Vol = nnz(MatC2) * prod(VoxDim);
            C2Vol_cm = C2Vol/1000;
                C2Vol_lts = C2Vol/100000;
        
        C3Vol = nnz(MatC3) * prod(VoxDim);
            C3Vol_cm = C3Vol/1000;
                C3Vol_lts = C3Vol/100000;
        
        [a,saida,c] = fileparts(filevol{i});

        if isequal(get(handles.IVR, 'Value'),1)   
            if isequal(get(handles.GM, 'Value'),1)
                nome = 'Grey matter';
                PartVol = C1Vol;
            end
            if isequal(get(handles.WM, 'Value'),1)  
                nome = 'White matter';
                PartVol = C2Vol;
            end
            if isequal(get(handles.CSF, 'Value'),1)
                nome = 'Cerebral Spinal Fluid';
                PartVol = C3Vol;
            end
            
            RatioVol = PartVol/TotVol;

            media2(i)=(RatioVol);

            Tex = sprintf('%s',saida);
            Tex1 = sprintf('The intracranial volume: %.3f cm³ ',TotVol/1000);
            Tex2 = sprintf('The %s volume is: %.3f cm³',nome,PartVol/1000);
            Tex3 = sprintf('The ratio between intracranial and %s volumes is: %.5f ',nome,RatioVol);

            fprintf(fide,'%s\r\n%s\r\n%s\r\n%s\r\n\r\n',Tex,Tex1,Tex2,Tex3);           
        else  
            if nOfSub>1
                media(i)=TotVol;
            end               
            Tex = sprintf('%s - Volume: %.2f mm³ %.3f cm³ or %.6f Litres',saida,TotVol,TotVol_cm,TotVol_lts);
            fprintf(fide,'%s\r\n\r\n',Tex);
        end
    end

    if nOfSub>1
        if get(handles.IVR, 'Value')==1  
            media2 = mean(media2);
            fprintf(fide,'\r\n');
            fprintf(fide,'The mean of ratios is: %.6f',media2);
        else
            media1 = mean(media);
            fprintf(fide,'The mean of volumes is: %.6f Litres',media1);
        end
    end
    fclose('all');
    set(handles.textresults,'Enable', 'on');
    set(handles.textshowed,'String', pathvol);
    set(handles.text12,'String', 'Done!!');
    set(handles.percent,'String', '');
              

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

function threefiles_Callback(hObject, eventdata, handles)
set(handles.Onefile, 'Value', 0)
set(handles.Run2,'Visible','off')

if get(handles.threefiles, 'Value') == 1
    set(handles.IVR,'Enable','on')
end

if get(handles.threefiles, 'Value') == 0
    set(handles.Onefile, 'Value', 1)
    set(handles.Run2,'Visible','on')
    set(handles.IVR,'Enable','off')
    set(handles.GM,'Enable','off')
    set(handles.WM,'Enable','off')
    set(handles.CSF,'Enable','off')
    set(handles.GM,'Enable','off')
    set(handles.GM,'Value',0)
    set(handles.WM,'Value',0)
    set(handles.CSF,'Value',0)
    set(handles.IVR,'Value',0)
end

function Onefile_Callback(hObject, eventdata, handles)
set(handles.threefiles, 'Value', 0)
set(handles.Run2,'Visible','on')
set(handles.IVR,'Enable','off')
set(handles.WM,'Enable','off')
set(handles.CSF,'Enable','off')
set(handles.GM,'Enable','off')
set(handles.GM,'Value',0)
set(handles.WM,'Value',0)
set(handles.CSF,'Value',0)
set(handles.IVR,'Value',0)

if get(handles.Onefile, 'Value') == 0
    set(handles.threefiles, 'Value', 1)
    set(handles.Run2,'Visible','off')
    set(handles.IVR,'Enable','on')
end

function Run2_Callback(hObject, eventdata, handles)
global filevol pathvol

    set(handles.text12,'String', 'Running.....Wait!');
    drawnow
    fid = fopen([pathvol filesep 'Volumetry_Results.txt'],'w+');
    AAA = size(filevol,2);
    
    for i = 1:AAA

        fprintf('Volume %d \n',i);
        
        Stru = nifti([pathvol filevol{i}]);
        matS = Stru.dat(:,:,:);

        x1 = Stru.hdr.pixdim(1,2);          
        y1 = Stru.hdr.pixdim(1,3);
        z1 = Stru.hdr.pixdim(1,4);
        PixVol = x1*y1*z1;
        
        matS(matS<str2num(get(handles.thrs,'String'))) = 0;
        matS(matS>=str2num(get(handles.thrs,'String'))) = 1;
        
        numVox = nnz(matS);
        mm3Vol = numVox*PixVol;
        CM = mm3Vol/1000;
        LTS = mm3Vol/1000000;

        Tex = sprintf('%s - Volume: %.2f mm³ %.5f cm³ or %.6f Litres / Num of Voxels: %d',filevol{1,i},mm3Vol,CM,LTS,numVox);
        fprintf(fid,'%s\r\n\r\n',Tex);
            
        if AAA>1
            media(i)=(LTS);
        end               
    end
    
    if AAA>1
        media1 = mean(media);
        fprintf(fid,'The mean volumes of this sample is: %.6f Litres',media1);
    end
    fclose(fid);
    set(handles.text12,'String', 'Done!');
    set(handles.textresults,'Enable', 'on');
    set(handles.textshowed,'String', pathvol);

% --- Executes on button press in IVR.
function IVR_Callback(hObject, eventdata, handles)

if get(handles.IVR,'Value')==1
    
    set(handles.GM,'Enable','on')
set(handles.WM,'Enable','on')
set(handles.CSF,'Enable','on')
set(handles.GM, 'Value',1)
else
    set(handles.GM,'Enable','off')
set(handles.WM,'Enable','off')
set(handles.CSF,'Enable','off')
set(handles.WM, 'Value',0)
set(handles.GM, 'Value',0)
set(handles.CSF, 'Value',0)
end


function WM_Callback(hObject, eventdata, handles)
set(handles.WM, 'Value',1)
set(handles.GM, 'Value',0)
set(handles.CSF, 'Value',0)

function GM_Callback(hObject, eventdata, handles)
set(handles.WM, 'Value',0)
set(handles.GM, 'Value',1)
set(handles.CSF, 'Value',0)


function CSF_Callback(hObject, eventdata, handles)
set(handles.WM, 'Value',0)
set(handles.GM, 'Value',0)
set(handles.CSF, 'Value',1)

function thrs_Callback(hObject, eventdata, handles)

function thrs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
