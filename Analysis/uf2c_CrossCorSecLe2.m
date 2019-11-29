function varargout = uf2c_CrossCorSecLe2(varargin)
% UF²C M-file for uf2c_CrossCorSecLe2.fig
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
                   'gui_OpeningFcn', @uf2c_CrossCorSecLe2_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_CrossCorSecLe2_OutputFcn, ...
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

function uf2c_CrossCorSecLe2_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_CrossCorSecLe2_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function push1_Callback(hObject, eventdata, handles)
global g1_VAR g2OK g1OK dirOK nodeOK

[fileG1,pathG1] = uigetfile({'*.mat','Matlab file'},...
        'Select the result file (All_Subjs-VAR.mat) for group 1','MultiSelect','off');
if ~isequal(fileG1,0)
    load([pathG1 fileG1])
    try
        g1_VAR = AllSubj3D;
    catch
        g1_VAR = TotalPairWise;
    end
    set(handles.text4,'String',sprintf('Data from %d subjects were added',size(g1_VAR,3)))
    g1OK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1)
        set(handles.runB,'Enable','on')
    end
end

function push2_Callback(hObject, eventdata, handles)
global g2_VAR g2OK g1OK dirOK nodeOK

[fileG2,pathG2] = uigetfile({'*.mat','Matlab file'},...
        'Select the result file (All_Subjs-VAR.mat) for group 2','MultiSelect','off');
if ~isequal(fileG2,0)
    load([pathG2 fileG2])
    try
        g2_VAR = AllSubj3D;
    catch
        g2_VAR = TotalPairWise;
    end
    set(handles.text5,'String',sprintf('Data from %d subjects were added',size(g2_VAR,3)))
    g2OK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1)
        set(handles.runB,'Enable','on')
    end
end

function COV_Callback(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK g1_VAR g2_VAR

if get(handles.COV,'Value')
    set(handles.G1cov,'Enable','on');
    set(handles.G2cov,'Enable','on');
else
    set(handles.G1cov,'Enable','off');
    set(handles.G2cov,'Enable','off');
end



function OutDir_Callback(hObject, eventdata, handles)
global Dir g2OK g1OK  dirOK nodeOK

Dir = uigetdir('','Define the output directory');
if ~isequal(Dir,0)
    Dir = [Dir filesep 'UF2C_SecondLevel'];
    set(handles.text10,'String',Dir)
    dirOK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1)
        set(handles.runB,'Enable','on')
    end
end

function runB_Callback(hObject, eventdata, handles)
global Dir COVDATA g1_VAR g2_VAR G1cov G2cov ROIdec
set(handles.text13,'String','Running.....')

mkdir(Dir)

if get(handles.COV,'Value')
    COVDATA = cat(3,G1cov,G2cov);
end

fprintf(1,'Calculating stats:     ');
for i = 1:size(g1_VAR,1)
    perc = round(100.*(i/size(g1_VAR,1)));
    if perc < 10
        fprintf(1,'\b\b%d%%',perc);
    else
        fprintf(1,'\b\b\b%d%%',perc);
    end
    for j = 1:size(g1_VAR,1)
        if i==j
            Pmap(i,j) = 1;
        else
%         [H,P,CI,tvalue] = ttest2(squeeze(g1_VAR(i,j,:)),squeeze(g2_VAR(i,j,:)));
            if get(handles.COV,'Value')
                [T,P,X1,X2] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                    [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],squeeze(COVDATA(i,j,:)),'group-group');
                 Pmap(i,j) = P(1);
%                  Tmap(i,j) = T(1);
            else
                [T,P,X1,X2] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                    [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],[],'group-group');
                Pmap(i,j) = P(1);
%                 Tmap(i,j) = T(1);
            end
        end
    end
end
fprintf('\n')
disp('-- Done,')
degF = size(g1_VAR,3) + size(g2_VAR,3);

% Uncorrected P Map
if get(handles.chec1,'Value')
    
    unCorPmap = Pmap;
    unCorPmap(unCorPmap>str2num(get(handles.edit1,'String'))) = 1;
    unCorPmap = -1.*unCorPmap;
    disp('Printing results - Uncorrected')
    fig1 = figure;
    imagesc(unCorPmap)
    title('Uncorrected P Map')
    colorbar
    caxis([-1.*(str2num(get(handles.edit1,'String'))) (max(max(unCorPmap))-(0.2*max(max(unCorPmap))))]) %reescala a volorbar para fugir dos extremos

    saveas(fig1,[Dir filesep 'Uncorrected_P_Map'],'png')
    saveas(fig1,[Dir filesep 'Uncorrected_P_Map'],'fig')
    save([Dir filesep 'unCorPmap-VAR'],'unCorPmap')
    close(fig1)

    % T map
    T_tmap = abs(unCorPmap);
    T_tmap(T_tmap==1)=NaN;
    T_tmap = abs(tinv(T_tmap,degF));
    Tvs = T_tmap(T_tmap>0);

    figA = figure;
    imagesc(T_tmap)
    title('Uncorrected T Map')
    colorbar
    caxis([(min(Tvs)-(0.2*min(Tvs))) max(Tvs)])

    saveas(figA,[Dir filesep filesep 'Uncorrected_T_Map'],'jpg')
    saveas(figA,[Dir filesep filesep 'Uncorrected_T_Map'],'fig')
    save([Dir filesep filesep 'Tmap-VAR'],'T_tmap')
    close(figA)

    if get(handles.checkbox9,'Value')
        summUnTmap = transpose(nansum(T_tmap));
        nx = 1;
        for rec = 1:size(ROIdec.DecVet,1)
            if ~isequal(ROIdec.DecVet(rec),0)
                NEWsumVetT(rec) = summUnTmap(nx);
                nx = nx+1;
            else
                NEWsumVetT(rec) = 0;
            end
        end

        PosROI3D = reshape(NEWsumVetT',ROIdec.MatS(1),ROIdec.MatS(2),ROIdec.MatS(3));
        PosROI3D_Stru = ROIdec.Stru;
        PosROI3D_Stru.dat.fname = [Dir filesep filesep 'Summed_T-Map.nii'];
        PosROI3D_Stru.dat.dtype = 'FLOAT32-LE';
        PosROI3D_Stru.dat(:,:,:) = PosROI3D;
        create(PosROI3D_Stru)
    end
    disp('-- Done,')
end

if get(handles.chec2,'Value')  %%%%% TCH corrected Matrix: Tukey, Ciminera and Heyse method
    disp('Printing results - Corrected')
    desiPvalue = str2num(get(handles.edit2,'String'));
    TCHcor = Pmap;
    TCHcor = tril(TCHcor,-1);
    AlphaC = 1-((1-desiPvalue)^(1/sqrt(size(g1_VAR,1)))); % Adjusted p-value threshold
    TCHcor(TCHcor>AlphaC)=NaN;
    TCHcor(TCHcor==0)=NaN;
    adj_p = TCHcor(TCHcor>0);
    TCHcor = -1.*TCHcor;
    
    if ~isequal(numel(unique(TCHcor)),1)
        fig2 = figure;
        imagesc(TCHcor)
        title('TCH corrected P Map')
        colorbar
        caxis([-1.*AlphaC (-1*(min(adj_p)-(0.2*min(adj_p))))])

        saveas(fig2,[Dir filesep filesep 'TCH_corrected_P_Map'],'jpg')
        saveas(fig2,[Dir filesep filesep 'TCH_corrected_P_Map'],'fig')
        save([Dir filesep filesep 'TCH_Pmap-VAR'],'TCHcor')
        close(fig2)

        % T map
        TCH_tmap = abs(TCHcor);
        TCH_tmap = abs(tinv(TCH_tmap,size(g1_VAR,1)));
        TCH_tmap(TCH_tmap==0)=NaN;
        adj_T = TCH_tmap(TCH_tmap>0);

        fig3 = figure;
        imagesc(TCH_tmap)
        title('TCH corrected T Map')
        colorbar
        caxis([(min(adj_T)-(0.2*min(adj_T))) max(adj_T)])

        saveas(fig3,[Dir filesep filesep 'TCH_corrected_T_Map'],'jpg')
        saveas(fig3,[Dir filesep filesep 'TCH_corrected_T_Map'],'fig')
        save([Dir filesep filesep 'TCH_Pmap-VAR'],'TCH_tmap')
        close(fig3)

         if get(handles.checkbox9,'Value')
            summCorTmap = nansum(TCH_tmap);
            nx = 1;
            for rec = 1:size(ROIdec.DecVet,1)
                if ~isequal(ROIdec.DecVet(rec),0)
                    NEWsumVetT2(rec) = summCorTmap(nx);
                    nx = nx+1;
                else
                    NEWsumVetT2(rec) = 0;
                end
            end

            PosROI3D = reshape(NEWsumVetT2',ROIdec.MatS(1),ROIdec.MatS(2),ROIdec.MatS(3));
            PosROI3D_Stru = ROIdec.Stru;
            PosROI3D_Stru.dat.fname = [Dir filesep filesep 'Summed_TCH_Correc_T-Map.nii'];
            PosROI3D_Stru.dat.dtype = 'FLOAT32-LE';
            PosROI3D_Stru.dat(:,:,:) = PosROI3D;
            create(PosROI3D_Stru)
         end
    else
        warndlg('No significant results with TCH correction','Ops!')
    end
    disp('-- Done,')
end
set(handles.text13,'String','Done!')
disp('All Done!')

function ROIfiles_Callback(hObject, eventdata, handles)
if get(handles.ROIfiles,'Value')
    set(handles.CoodList,'Value',0)
    set(handles.addNODES,'String','Add ROIs images')
else
    set(handles.CoodList,'Value',1)
    set(handles.addNODES,'String','Add Coordinate list')
end

function chec1_Callback(hObject, eventdata, handles)
function chec2_Callback(hObject, eventdata, handles)
function CorMats_Callback(hObject, eventdata, handles)

function Grapho2D_Callback(hObject, eventdata, handles)
global g2OK g1OK  dirOK nodeOK

if isequal(get(handles.Grapho2D,'Value'),0)
    if logical(g2OK) && logical(g1OK) && logical(dirOK)
        if isequal(get(handles.Graphos3D,'Value'),0)
            set(handles.runB,'Enable','on')
        end
    end
else
    if logical(g2OK) && logical(g1OK) && logical(dirOK) && logical(nodeOK)
        set(handles.runB,'Enable','on')
    else
        set(handles.runB,'Enable','off')
    end
end

function Graphos3D_Callback(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK

if isequal(get(handles.Graphos3D,'Value'),0)
    if logical(g2OK) && logical(g1OK) && logical(dirOK)
        if isequal(get(handles.Grapho2D,'Value'),0)
            set(handles.runB,'Enable','on')
        end
    end
else
    if logical(g2OK) && logical(g1OK) && logical(dirOK) && logical(nodeOK)
        set(handles.runB,'Enable','on')
    else
        set(handles.runB,'Enable','off')
    end
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


function figure1_CreateFcn(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK
g2OK=0; g1OK=0; dirOK=0; nodeOK=0;


function G1cov_Callback(hObject, eventdata, handles)
global G1cov

[G1cov,pathG1cov] = uigetfile({'*.mat','Matlab file'},...
        'Select a .mat file for covariate','MultiSelect','off');
load([pathG1cov G1cov]);
G1cov = All_Subjs_DiffGM;
clear All_Subjs_DiffGM

function G2cov_Callback(hObject, eventdata, handles)
global G2cov

[G2cov,pathG2cov] = uigetfile({'*.mat','Matlab file'},...
        'Select a .mat file for covariate','MultiSelect','off');
load([pathG2cov G2cov]);
G2cov = All_Subjs_DiffGM;
clear All_Subjs_DiffGM

function checkbox9_Callback(hObject, eventdata, handles)

if get(handles.checkbox9,'Value')
    set(handles.pushbutton10,'Enable','on')
else
    set(handles.pushbutton10,'Enable','off')
end

function pushbutton10_Callback(hObject, eventdata, handles)
global ROIdec

[ROIdec,pathROIdec] = uigetfile({'*.mat','Matlab file'},...
        'Select the ROIdecoderResh.mat file of the used ROI','MultiSelect','off');
load([pathROIdec ROIdec]);
ROIdec = ROI_Decoder;
