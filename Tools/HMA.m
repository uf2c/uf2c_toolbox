function varargout = HMA(varargin)
% UF²C M-file for HMA.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2015
%
% Copyright (c) 2015, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @HMA_OpeningFcn, ...
                   'gui_OutputFcn',  @HMA_OutputFcn, ...
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
function HMA_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);
function varargout = HMA_OutputFcn(hObject, eventdata, handles) 


function AddRP_Callback(hObject, eventdata, handles)
global filetxt

filetxt = uipickfiles('Output','cell','Prompt','Add all subjects files','FilterSpec','rp_*.txt');

set(handles.AddRPtxt,'String',[num2str(size(filetxt,2)) ' rp file(s) added'])

function runB_Callback(hObject, eventdata, handles)
global filetxt OutDir
set(handles.text11,'String','Running...')
drawnow

% dateNow = clock;
% Fname = (['1-RESULTs_Mov_Check_' num2str(dateNow(3)) '_' num2str(dateNow(2)) '_' num2str(dateNow(1)) '--' num2str(dateNow(4)) '_' num2str(dateNow(5)) '.txt']);

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

ExOK = 0;
[sxt,dxt,fxt] = fileparts(filetxt{1,1});

if isempty(OutDir)
    OutDir = [sxt filesep];
else
    OutDir = [OutDir filesep];
end

fideG = fopen([OutDir,get(handles.Resname,'String'),'.txt'],'w+'); % CHARGE OUTPUT LOG FILE
fprintf(fideG,'                                    Head Motion Analysis\r\n\r\n');
fprintf(fideG,'Number of subjects included: \t%d \r\n\r\n', size(filetxt,2));
fprintf(fideG,...
    'Case \t Max_X \t Max_Y \t Max_Z \t Max_Amplitude_T \t Max_Tranlation \t Sum_Translation \t Max_Row \t Max_Pitch \t Max_Yaw \t Max_Amplitude_Rot \t Max_Rotation \t Sum_Rotation \t Framewise Displacement (FD)\r\n');

asu = 1;
for i = 1:size(filetxt,2)
    
    [sx,dx,fx] = fileparts(filetxt{1,i});
    vetore = importdata(filetxt{1,i});
    vetore(:,4:6) = (180/pi).*vetore(:,4:6);
    
    if get(handles.checkbox5,'Value')
        [Mot_fig,HDvV,HTvV] = uf2c_plot_motion(filetxt{1,i},'off');
        imgRR = getframe(Mot_fig);
        imwrite(imgRR.cdata,[sx,filesep,dx,'_RP_Plot.png']);
        saveas(Mot_fig,[sx,filesep,dx,'_RP_Plot'],'fig');
        close(Mot_fig)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BRAMILA's Framewise Displacemente (FD)
    % Code from: Brain and Mind Lab at Aalto University
    % Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also 
    % Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg.motionparam = filetxt{1,i};
    cfg.prepro_suite = 'spm';
    cfg.radius = 50;

    [FDts,rms] = bramila_framewiseDisplacement(cfg);

%     save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement.mat'],'FDts')
    if get(handles.checkbox4,'Value')
        figuFD = figure('Visible','off');
        set(figuFD,'Name','Framewise Displacemente (FD) TS',...
            'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
            'Color',[1 0.94 0.86]);
        plot(FDts);
        hold on
        plot(ones(1,numel(FDts)).*str2num(get(handles.edit4,'String')));
        ylabel('FD (mm)')
        xlabel('Time Points')
        title('Average Framewise Displacement TS','FontSize', 14);
        drawnow    
        imgRR = getframe(figuFD);
        imwrite(imgRR.cdata, [sx,filesep,dx,'_FD_Plot.png']);
        saveas(figuFD,[sx,filesep,dx,'_FD_Plot'],'fig')
        close(figuFD)
    end

    vetX = vetore(:,1);
    vetY = vetore(:,2);
    vetZ = vetore(:,3);
    vetR = vetore(:,4);
    vetP = vetore(:,5);
    vetYa = vetore(:,6);
    
    MaxX = max(abs(vetX));
    MaxY = max(abs(vetY));
    MaxZ = max(abs(vetZ));
    sumT = sum(sum(abs(vetore(:,1:3))));
    MaxT = max([MaxX,MaxY,MaxZ]);
    
    MaxR = max(abs(vetR));
    MaxP = max(abs(vetP));
    MaxYa = max(abs(vetYa));
    sumRot = sum(sum(abs(vetore(:,4:6))));
    MaxRot = max([MaxR,MaxP,MaxYa]);

    for j = 1:size(vetore,1)-1
        vetX(j) = abs(vetX(j)-vetX(j+1));
        vetY(j) = abs(vetY(j)-vetY(j+1));
        vetZ(j) = abs(vetZ(j)-vetZ(j+1));
        vetR(j) = abs(vetR(j)-vetR(j+1));
        vetP(j) = abs(vetP(j)-vetP(j+1));
        vetYa(j) = abs(vetYa(j)-vetYa(j+1));
    end
    
    vetX = vetX(1:end-1);
    vetY = vetY(1:end-1);
    vetZ = vetZ(1:end-1);
    MAXaMP_t = max([max(vetX),max(vetY),max(vetZ)]);
    
    vetR = vetR(1:end-1);
    vetP = vetP(1:end-1);
    vetYa = vetYa(1:end-1);
    MAXaMP_r = max([max(vetR),max(vetP),max(vetYa)]);
    
    fprintf(fideG,'%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t\r\n',...
                   dx,MaxX,MaxY,MaxZ,MAXaMP_t,MaxT,sumT,MaxR,MaxP,MaxYa,MAXaMP_r,MaxRot,sumRot,mean(FDts));
    
    if get(handles.checkbox1,'Value') && get(handles.checkbox2,'Value') && get(handles.checkbox3,'Value') 
        if  MaxT>str2num(get(handles.Tthresh,'String')) || MaxRot>str2num(get(handles.Rthresh,'String')) || mean(FDts)>str2num(get(handles.edit4,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if get(handles.checkbox1,'Value') && get(handles.checkbox2,'Value') && isequal(get(handles.checkbox3,'Value'),0)
        if  MaxT>str2num(get(handles.Tthresh,'String')) || MaxRot>str2num(get(handles.Rthresh,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if get(handles.checkbox1,'Value') && isequal(get(handles.checkbox2,'Value'),0) && get(handles.checkbox3,'Value') 
        if  MaxT>str2num(get(handles.Tthresh,'String')) || mean(FDts)>str2num(get(handles.edit4,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if isequal(get(handles.checkbox1,'Value'),0) && get(handles.checkbox2,'Value') && get(handles.checkbox3,'Value') 
        if  MaxRot>str2num(get(handles.Rthresh,'String')) || mean(FDts)>str2num(get(handles.edit4,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if get(handles.checkbox1,'Value') && isequal(get(handles.checkbox2,'Value'),0) && isequal(get(handles.checkbox3,'Value'),0)
        if  MaxT>str2num(get(handles.Tthresh,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if isequal(get(handles.checkbox1,'Value'),0) && get(handles.checkbox2,'Value') && isequal(get(handles.checkbox3,'Value'),0)
        if  MaxT>str2num(get(handles.Rthresh,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end
    
    if isequal(get(handles.checkbox1,'Value'),0) && isequal(get(handles.checkbox2,'Value'),0) && get(handles.checkbox3,'Value')
        if  mean(FDts)>str2num(get(handles.edit4,'String'))
            Exc(asu) = i;
            asu = asu+1;
            ExOK = 1;
        end
    end

end

if ExOK
    fprintf(fideG,'\r\n');
    fprintf(fideG,'Excluded due excessive movement (relative to defined thresholds)\r\n');
    for t = 1:asu-1
        [de name ext] = fileparts(filetxt{1,Exc(t)});
        fprintf(fideG,'%s\r\n',name);
    end
end
fclose('all');
disp('Done!')
set(handles.text11,'String','Done!')

function Tthresh_Callback(hObject, eventdata, handles)

function edit4_Callback(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)
if get(handles.checkbox1,'Value')
    set(handles.text5,'Enable','on')
    set(handles.Tthresh,'Enable','on')
else
    set(handles.text5,'Enable','off')
    set(handles.Tthresh,'Enable','off')
end

function checkbox2_Callback(hObject, eventdata, handles)
if get(handles.checkbox2,'Value')
    set(handles.text6,'Enable','on')
    set(handles.Rthresh,'Enable','on')
else
    set(handles.text6,'Enable','off')
    set(handles.Rthresh,'Enable','off')
end

function checkbox3_Callback(hObject, eventdata, handles)
if get(handles.checkbox3,'Value')
    set(handles.text8,'Enable','on')
    set(handles.edit4,'Enable','on')
    set(handles.checkbox4,'Enable','on')
else
    set(handles.text8,'Enable','off')
    set(handles.edit4,'Enable','off')
    set(handles.checkbox4,'Enable','off')
    set(handles.checkbox4,'Value',0)
end

function Rthresh_Callback(hObject, eventdata, handles)

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)

function OutDir_Callback(hObject, eventdata, handles)
global OutDir

OutDir = uigetdir('','Define the output directory');

if isequal(OutDir,0)
    clear OutDir;
else
    set(handles.text9,'String',OutDir)
    drawnow
end

function Resname_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
