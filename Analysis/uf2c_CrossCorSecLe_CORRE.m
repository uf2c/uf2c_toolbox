function varargout = uf2c_CrossCorSecLe_CORRE(varargin)
% UF²C M-file for uf2c_CrossCorSecLe_CORRE.fig
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
                   'gui_OpeningFcn', @uf2c_CrossCorSecLe_CORRE_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_CrossCorSecLe_CORRE_OutputFcn, ...
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

function uf2c_CrossCorSecLe_CORRE_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_CrossCorSecLe_CORRE_OutputFcn(hObject, eventdata, handles) 
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
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function COV_Callback(hObject, eventdata, handles)
global COVDATA nOfLines nOfCovs g1OK dirOK nodeOK g1_VAR CovHeaders

if get(handles.COV,'Value')
    [COVlist,COVpath] = uigetfile({'*.txt','Text files'},'Select a text file with all covariates','MultiSelect','off');
    if ~isequal(COVlist,0)
        COVDATA = importdata([COVpath COVlist]);
        if isstruct(COVDATA)
            nOfLines = size(COVDATA.data,1);
            nOfCovs = size(COVDATA.data,2);
            CovHeaders = COVDATA.colheaders;
            COVDATA = COVDATA.data;
        else
            nOfLines = size(COVDATA,1);
            nOfCovs = size(COVDATA,2);
            CovHeaders = strsplit(num2str(1:nOfCovs));
        end
        if isequal(g1OK,1)
            if isequal(size(g1_VAR,3),nOfLines)
                set(handles.text9,'String',sprintf('%d Correlation variable(s) added!',nOfCovs))
            else
                warndlg('The number of inputs in each vector column should to match with the Group n of subjects','Ops!')
                set(handles.COV,'Value',0)
                set(handles.text9,'String','')
            end
        else
            set(handles.text9,'String',sprintf('%d Correlation vector(s) added!',nOfCovs))
        end
        if isequal(g1OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    else
        set(hanldes.COV,'Value',0)
        set(hanldes.text9,'String','')
    end
else
    set(handles.text9,'String','')
end

function OutDir_Callback(hObject, eventdata, handles)
global Dir g1OK  dirOK nodeOK

Dir = uigetdir('','Define the output directory');
if ~isequal(Dir,0)
    Dir = [Dir filesep 'UF2C_SecondLevel_Correlations'];
    set(handles.text10,'String',Dir)
    dirOK = 1;
    if isequal(g1OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function addNODES_Callback(hObject, eventdata, handles)
global CoordMS seedList nodeOK g1OK dirOK nEle nOfNets

if get(handles.CoodList,'Value')
    set(handles.text12,'String','Wait.......')
    [fileCor,pathCor] = uigetfile({'*.txt','Text file'},...
            'Select the cordinates list file','MultiSelect','off');
        if ~isequal(fileCor,0)
            seedListT = importdata([pathCor fileCor]);
            nOfNets = str2num(seedListT{end,1}(1:2));
            xs=1;
            nE = 1;
            while xs <= size(seedListT,1)
                xStr = seedListT{xs,1}(1:3);
                net1 = strncmp(xStr,seedListT,3);
                nEle(nE) = sum(net1);
                xs = xs + sum(net1);
                nE = nE+1;
            end
            for i = 1:size(seedListT,1)
                seedList(i,:) = str2num(seedListT{i}(4:end));
            end
            RootDir = which('uf2c');
            RootDir = ([RootDir(1:end-6) 'Image_Examples' filesep]);
            for jh = 1:size(seedList,1)
                Vhdr = spm_vol([RootDir,'FiltRegrSW_Example.nii']);
                TcM = Vhdr.mat;
                EPIcoord = [seedList(jh,1) seedList(jh,2) seedList(jh,3) 1]*(inv(TcM))';
                EPIcoord = round(EPIcoord);
                CoordMS(jh,:) = EPIcoord(1,:);
            end
            set(handles.text12,'String',sprintf('%d coordinates of %d networks were added!',size(CoordMS,1),nOfNets))
            nodeOK = 1;
            if isequal(g1OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
                set(handles.runB,'Enable','on')
            end
        end
else
    set(handles.text12,'String','Wait.......')
    [fileROI,pathROI] = uigetfile({'*.nii','NIfTI files'},'Select 3D ROI masks',...
        'MultiSelect','on');
    if ~isequal(fileROI,0)
        if ~iscell(fileROI)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
            fileROI = {fileROI};
        end
        fileROI = sort(fileROI);
        xs=1;
        nE = 1;
        loopSize = size(fileROI,2);
        while xs <= loopSize
           xStr = fileROI{1,xs}(1:3);
           net1 = strncmp(xStr,fileROI,3);
           nEle(nE) = sum(net1);
           xs = xs + sum(net1);
           nE = nE+1;
        end
        nOfNets = size(nEle,2);

        for jh = 1:loopSize
            struRoi = nifti([pathROI fileROI{jh}]);
            ROIMAt = struRoi.dat(:,:,:);
            ROIMAt(ROIMAt>0)=1;
            roi1 = ROIMAt;
%             [iCo,jCo,kCo]= ind2sub(size(roi1), find(roi1>0));
%             coords = [iCo,jCo,kCo];
            ncentr = 1;
            for yt = 1:size(roi1,3)
                roiFind = roi1(:,:,yt);
                roiFind = fliplr(roiFind);
                roiFind = permute(roiFind,[2,1]);
                roiFind = fliplr(roiFind);
                roiFind = flipud(roiFind);
                if ~isequal(sum(sum(roiFind)),0)
                    roiFind2 = logical(roiFind);
                    s = regionprops(roiFind2,'centroid');
                    if size(s,1)>1
                        for ytr = 1:size(s,1)
                            centNs(ytr,:) = s(ytr).Centroid;
                        end
                        sFnl(ncentr,1:2) = mean(centNs,1);
                        sFnl(ncentr,3) = yt;  
                    else
                        sFnl(ncentr,1:2) = s.Centroid(1,:);
                        sFnl(ncentr,3) = yt;  
                    end
                    ncentr = ncentr+1;
                    clear s
                end
            end
            FinalCentrCoord(jh,:) = round(median(sFnl,1));
            clear sFnl
        end
        nodeOK = 1;
        CoordMS = FinalCentrCoord;
        seedList = transpose(fileROI);
        set(handles.text12,'String',sprintf('%d coordinates of %d netorworks were added!',size(CoordMS,1),nOfNets))
        if isequal(g1OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    end
end
    
function runB_Callback(hObject, eventdata, handles)
global CoordMS seedList Dir COVDATA g1_VAR nOfLines nEle nOfNets CovHeaders
set(handles.text13,'String','Running.....')
%%%%%% Checking...
if ~isequal(size(CoordMS,1),size(g1_VAR,1))
    if logical(get(handles.Grapho2D,'Value')) || logical(get(handles.Graphos3D,'Value'))
        warndlg(sprintf('The number of nodes added (%d) should to match with the used on first level (%d)',size(CoordMS,1),size(g1_VAR,1)),'Ops! Process Aborted')
    return
    end
end
if get(handles.COV,'Value')
    if ~isequal(size(g1_VAR,3),nOfLines)
        warndlg(sprintf('The number of data of each covariate (%d) should to match with the sum of Group 1 and 2 subjects',nOfLines,(size(g1_VAR,3))),'Ops! Process Aborted')
        return
    end
end
%%%%%%%%%%%%%%%%
mkdir(Dir)

NetLabels = strcat('n',strsplit(num2str(1:nOfNets)));
DirF = Dir;
    g1_VAR2 = g1_VAR;
    fact = 1;
    TarNetStr = [''];
    TarNetC = [];

for x = 1:size(COVDATA,2)
    fprintf('Processing the variable ''%s''',CovHeaders{x});
    
    for i = 1:size(g1_VAR,1)
        for j = 1:size(g1_VAR,1)
            if i==j
                Pmap(i,j) = 1;
                Rmap(i,j) = 0;
            else
                [R,P] = corrcoef(squeeze(g1_VAR(i,j,:)),COVDATA(:,x));
                 Pmap(i,j) = P(2);
                 Rmap(i,j) = R(2);
            end
        end
    end
    Rmap(isnan(Rmap)) = 0;
    Pmap(isnan(Pmap)) = 1;
    ScreSize = get(0,'screensize');
    ScreSize = ScreSize(3:end);

    fig00 = figure;
    imagesc(Rmap)
    set(fig00,'Name','Correlation Matrix','Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
            'Color',[1 .949 .867]);
    title('Correlation Matrix','FontSize',14);
    xlabel('Seeds'  ,'FontSize',14);
    ylabel('Seeds','FontSize',14);
    set(gca, 'color', [0 0 0])
    try  set(gca, 'gridcolor', [1 1 1]); end
    try set(gca, 'gridalpha', 1); end
    daspect([1 1 1])
    colorbar
    caxis([min(min(Rmap)) max(max(Rmap))]) %reescala a colorbar para fugir dos extremos

    drawnow
    imgRR = getframe(fig00);
    imwrite(imgRR.cdata, [Dir filesep 'Correlation_Map_' CovHeaders{x} '.png']);
    clear imgRR
    
    saveas(fig00,[Dir filesep 'Correlation_Map_' CovHeaders{x}],'fig')
    save([Dir filesep 'Rmap'],'Rmap')
    close(fig00)

    stIndx = 1;
    for hg = 1:size(nEle,2)
        str_nEle(hg) = stIndx;
        stIndx = stIndx+nEle(hg);
    end
    g1_mean = mean(g1_VAR,3);


    
% Uncorrected P Map
if get(handles.chec1,'Value')
    
    unCorPmap = Pmap;
    unCorPmap(unCorPmap>str2num(get(handles.edit1,'String'))) = 1;
    unCorPmap = -1.*unCorPmap;

    if ~isequal(numel(unique(unCorPmap)),1)
        if get(handles.CorMats,'Value')   %%%%% Correlation Matrix
            fig1 = figure;
            imagesc(unCorPmap)

            set(fig1,'Name','Uncorrected Correlation Matrix','Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                    'Color',[1 .949 .867]);
            title(['Uncorrected Correlation Matrix (' CovHeaders{x} ')'],'FontSize',14);
            xlabel('Seeds'  ,'FontSize',14);
            ylabel('Seeds','FontSize',14);
            set(gca, 'color', [0 0 0])
            try  set(gca, 'gridcolor', [1 1 1]); end
            try set(gca, 'gridalpha', 1); end
            daspect([1 1 1])
            colorbar
            caxis([-1.*(str2num(get(handles.edit1,'String'))) (max(max(unCorPmap))-(0.2*max(max(unCorPmap))))]) %reescala a volorbar para fugir dos extremos
            
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_R_Map_' CovHeaders{x} '.png']);
            clear imgRR
            saveas(fig1,[Dir filesep 'Uncorrected_R_Map_' CovHeaders{x}],'fig')
            save([Dir filesep 'unCorRmap-VAR'],'unCorPmap')
            close(fig1)
        end
        matrixR = unCorPmap;
        matrixR(matrixR==-1) = 0;
        matrixR = abs(matrixR);

        % Creating Output text file
        matrix_bin = matrixR>0;

        fidRes = fopen([DirF,filesep,'Uncorr_Results_Description_' CovHeaders{x} '.txt'],'w+');
        
        for r = 1:size(g1_VAR,1)
            Ridx = find(matrix_bin(r,:));
            if sum(Ridx)>0
                tester1 = r<=(str_nEle+nEle)-1;
                networkNumb1 = find(tester1==1,1);
                try
                    RnUmb1 = r - sum(nEle(1:networkNumb1-1));
                catch
                    RnUmb1 = r;
                end
                fprintf(fidRes,'ROI %d: r%d%s: \r\n',r,RnUmb1,NetLabels{networkNumb1});
                for ixR = 1:numel(Ridx)
                     if Rmap(r,Ridx(ixR))>0
                         strDir = 'Direct Correlation';
                     else
                         strDir = 'Inverse Correlation';
                     end
                    tester2 = Ridx(ixR)<=(str_nEle+nEle)-1;
                    networkNumb2 = find(tester2==1,1);
                    try
                        RnUmb2 = Ridx(ixR) - sum(nEle(1:networkNumb2-1));
                    catch
                        RnUmb2 = Ridx(ixR);
                    end
                    fprintf(fidRes,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),RnUmb2,NetLabels{networkNumb2},strDir);
                end
            end
        end
        fclose(fidRes);
        
        limiteG = (size(g1_VAR,2)/2).*1.1;

        RAIO = limiteG/2.2;
        ang=0:0.01:2*pi; 
        xp=RAIO*cos(ang);
        yp=RAIO*sin(ang);
        xpt=(limiteG/2.15)*cos(ang);
        ypt=(limiteG/2.15)*sin(ang);

        AnoRa = size(xp,2)/size(matrixR,2);
        xp2 = xp(1:AnoRa:end);
        yp2 = yp(1:AnoRa:end);

        xp2t = xpt(1:AnoRa:end);
        yp2t = ypt(1:AnoRa:end);
        xp2t = xp2t+(limiteG/2);
        yp2t = yp2t+(limiteG/2);
        xp2t = xp2t(end:-1:1);
        yp2t = yp2t(end:-1:1);

        xp2 = xp2+(limiteG/2);
        yp2 = yp2+(limiteG/2);
        yp2 = yp2(end:-1:1);
        xp2 = xp2(end:-1:1);
        
        minH = (pdist([xp2(1),yp2(1);xp2(2),yp2(2)],'euclidean'));
        
        stIndx = 1;

         for hg = 1:size(nEle,2)
             str_nEle(hg) = stIndx;
             stIndx = stIndx+nEle(hg);
         end

         cc=hsv(nOfNets);
         ccV1 = cc(1:2:end,:);
         ccV2 = cc(2:2:end,:);
         cc = [ccV1;ccV2];

        LoopAte = size(matrixR,1);
        TranspMin = 0;
        Diam = max(xp)-min(xp);
        
        circ = figure('Visible','on');
        set(circ,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                  'Color',[0 0 0]);
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Circular Connectome - DRAWING........');

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ')']);
        drawnow

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);
            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end
            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                    h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                    R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                    hold on
                    NewC = circ_cent(P1,P2,R);

                    Marg = (limiteG-2*RAIO)/2;
                    DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                    DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                    if DisC1>DisC2
                        centerIdx = 1;
                    else
                        centerIdx = 2;
                    end

                    AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                    AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                    if sum(abs(AngRad))>pi && any(AngRad<0)
                        idxCic = find(AngRad<0);
                        AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                    end

                    [AngRad,ksord] = sort(AngRad);

                    if isequal(ksord,[1 2])
                        TranOrd = 'descend';
                    else
                        TranOrd = 'ascend';
                    end

                    ang2 = AngRad(1):0.03:AngRad(2);
                    ang2 = [ang2,AngRad(2)];

                    xunit = R * cos(ang2) + NewC(centerIdx,1);
                    yunit = R * sin(ang2) + NewC(centerIdx,2);

                    TrnaspSca = TranspMin:(1-TranspMin)/size(xunit,2):1;
                    TrnaspSca = sort(TrnaspSca,TranOrd);
                    figure(circ)
                    for hg = 1:size(xunit,2)-1
                        h = plot(xunit(1,hg:hg+1),yunit(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                        h.Color(4) = TrnaspSca(hg);
                    end
                end
            end
        end
        set(mTextBox1,'String',['Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circ, 'InvertHardCopy', 'off');
        imgRR = getframe(circ);
        imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Circ_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circ)

        circTS = figure('Visible','on');
        set(circTS,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                  'Color',[0 0 0],'Units','Normalized');
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        pos1 = [0.15 0.15 0.7 0.7];
        subplot('Position',pos1,'Color',[0 0 0]);
        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);
        axis off

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Circular Connectome - DRAWING........');

        TvaCrit = tinv(1-(str2num(get(handles.edit1,'String'))*.5),(size(g1_VAR2,3)-1));

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.3  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ', |Tscore|>' num2str(TvaCrit) ')']);
        drawnow

        ccJet = jet(64);
        ccGre = gray(64);
        for xcv = 1:size(g1_VAR2,1)
            for ycv = 1:size(g1_VAR2,2)
                Tva(xcv,ycv) = tinv(1-(matrixR(xcv,ycv)*.5),(size(g1_VAR2,3)-1));
            end
        end
        Tva(Tva==inf) = 0;
        TvaMax = max(max(Tva));

        cbm1 = 0;
        cbm2 = 0;
        cbm3 = 0;

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);

            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end

            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    if Rmap(i-(fact-1),vet(j))>0
                        colorVlidx = (round((Tva(i-(fact-1),vet(j))/TvaMax)*32))+32;
                        colorVl = ccJet(colorVlidx,:);
                    else
                        colorVlidx = round(((-Tva(i-(fact-1),vet(j))/TvaMax)*32))+33;
                        colorVl = ccJet(colorVlidx,:);
                    end

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                    h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                    R = ((ARCalpha^2) + (4*(h^2)))/(8*h);

                    hold on
                    NewC = circ_cent(P1,P2,R);
                    Marg = (limiteG-2*RAIO)/2;
                    DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                    DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                    if DisC1>DisC2
                        centerIdx = 1;
                    else
                        centerIdx = 2;
                    end

                    AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                    AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                    if sum(abs(AngRad))>pi && any(AngRad<0)
                        idxCic = find(AngRad<0);
                        AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                    end

                    [AngRad,ksord] = sort(AngRad);

                    if isequal(ksord,[1 2])
                        TranOrd = 'descend';
                    else
                        TranOrd = 'ascend';
                    end

                    ang2 = AngRad(1):0.03:AngRad(2);
                    ang2 = [ang2,AngRad(2)];

                    xunit = R * cos(ang2) + NewC(centerIdx,1);
                    yunit = R * sin(ang2) + NewC(centerIdx,2);
                    figure(circTS)
                    h1 = plot(xunit,yunit,'Color',colorVl,'LineWidth',3);
                    h1.Color(4) = 0.4;
                end
            end
        end

        clear A_cell Labx1
        subplot('Position',[0.04 0.01 0.9 0.015])
        imagesc(-TvaMax:(TvaMax/32):TvaMax);
        colormap(gca,'jet');
        set(gca, 'XAxisLocation', 'top','YColor',[0 0 0],'XColor',[1 1 1])
        xticks(1:64/16:64);
        Labx1 = [(round(10*(-TvaMax:(TvaMax/8):TvaMax)))/10 TvaMax];
        A_cell = split(cellstr(num2str(Labx1)));
        xticklabels(A_cell)
        yticklabels('')
        title(gca,'T-Score Direct Corres. (hot colors) Inverse Corres. (cool colors)','Color','w')

        set(mTextBox1,'String',['Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circTS, 'InvertHardCopy', 'off');
        imgRR = getframe(circTS);
        imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_T-Score_Circ_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circTS)
        
        circ2 = figure('Visible','on');
            set(circ2,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0]);
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.3  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Organic Circular Connectome - DRAWING........');

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ')']);
        drawnow

        daspect([1 1 1])
        RAIO2 = RAIO*0.15;        
        xp3=RAIO2*cos(ang);
        yp3=RAIO2*sin(ang);

        xp3 = xp3(1:AnoRa:end);
        yp3 = yp3(1:AnoRa:end);
        xp3 = flip(xp3+(limiteG/2));
        yp3 = flip(yp3+(limiteG/2));

        RAIO3 = RAIO*0.85;        
        xp4=RAIO3*cos(ang);
        yp4=RAIO3*sin(ang);

        idxNR = ceil((str_nEle+(nEle./2))-1);
        xp4 = flip(xp4(1:AnoRa:end));
        yp4 = flip(yp4(1:AnoRa:end));
        xp4 = xp4(idxNR);
        yp4 = yp4(idxNR);

        xp4 = xp4+(limiteG/2);
        yp4 = yp4+(limiteG/2);

        RAIO4 = RAIO*0.95;        
        xp5=RAIO4*cos(ang);
        yp5=RAIO4*sin(ang);

        xp5 = xp5((AnoRa/2):AnoRa:end);
        yp5 = yp5((AnoRa/2):AnoRa:end);
        xp5 = flip(xp5+(limiteG/2));
        yp5 = flip(yp5+(limiteG/2));

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);
            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end
            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    testerj = vet(j)<=(str_nEle+nEle)-1;
                    networkNumbj = find(testerj==1,1);
                    colorVj = cc(networkNumbj,:);
                    xxF = [];
                    yyF = [];
                    if networkNumb~= networkNumbj
                        opt1 =  round(i-(i-vet(j))/2);
                        opt2 = round(opt1+(LoopAte/2));
                        if opt2<=LoopAte
                            dopt1 = pdist([xp3(opt1),yp3(opt1);P1(1),P1(2)]);
                            dopt2 = pdist([xp3(opt2),yp3(opt2);P1(1),P1(2)]);
                            if dopt1>dopt2
                                dopt = opt2;
                            else
                                dopt = opt1;
                            end
                        else
                            dopt = opt1;
                        end

                        [xx,yy]=fillline([P1(1),P1(2)],[xp4(networkNumb),yp4(networkNumb)],15);
                        [xx2,yy2]=fillline([xp4(networkNumb),yp4(networkNumb)],[xp3(dopt),yp3(dopt)],20);
                        [xx3,yy3]=fillline([xp3(dopt),yp3(dopt)],[xp4(networkNumbj),yp4(networkNumbj)],20);
                        [xx4,yy4]=fillline([xp4(networkNumbj),yp4(networkNumbj)],[P2(1),P2(2)],15);

                        xxF = [xx,xx2,xx3,xx4];
                        yyF = [yy,yy2,yy3,yy4];
                        xxF = [xxF(1);smoothLNI(xxF(2:end-1),30);xxF(end)]';
                        yyF = [yyF(1);smoothLNI(yyF(2:end-1),30);yyF(end)]';
                        TranOrd = 'descend';

                        TrnaspSca = [zeros(1,size(xxF,2)*0.3)+TranspMin,TranspMin:(1-TranspMin)/(size(xxF,2)*0.7):1];
                        TrnaspSca = sort(TrnaspSca,TranOrd);
                        figure(circ2)
                        if networkNumb~=networkNumbj
                            for hg = 1:(size(xxF,2)-1)*0.7
                                h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                h.Color(4) = TrnaspSca(hg);
                            end
                        else
                            if i<vet(j)
                                h = plot(xxF,yyF,'Color',colorV,'LineWidth',1);
                            end
                        end

                    else

                        ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                        h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                        R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                        NewC = circ_cent(P1,P2,R);
                        Marg = (limiteG-2*RAIO)/2;
                        DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                        DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                        if DisC1>DisC2
                            centerIdx = 1;
                        else
                            centerIdx = 2;
                        end

                        AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                        AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                        if sum(abs(AngRad))>pi && any(AngRad<0)
                            idxCic = find(AngRad<0);
                            AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                        end

                        ang2 = AngRad(1):0.03:AngRad(2);
                        ang2 = [ang2,AngRad(2)];

                        xunit = R * cos(ang2) + NewC(centerIdx,1);
                        yunit = R * sin(ang2) + NewC(centerIdx,2);

                        h = plot(xunit,yunit,'Color',colorV.*0.5,'LineWidth',2);
                    end

                    if vet(j)<10
                        strf = ['0' num2str(vet(j))];
                    else
                        strf = num2str(vet(j));
                    end
                    bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                        colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                        'VerticalAlignment','middle');

                    Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                end
            end
        end

        set(mTextBox1,'String',['Organic Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circ2, 'InvertHardCopy', 'off');
    %     saveas(circ,[DirF filesep 'Uncorrected_Circ_Map'],'fig')
        imgRR = getframe(circ2);
        imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_OrganicCirc_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circ2)
            
        if get(handles.Graphos3D,'Value')  %%%%% 3D Graphos
            GraphoFig = figure;

            set(GraphoFig,'Name','3DConnectome',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                            'Color',[1 .949 .867]);
                        
            title('Uncorrected 3D Connectome','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
%             scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0.8 0.8 0.8]); end
            
            try set(gca, 'gridalpha', 1); end
            
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            view([-1 -1 1])
            hold on

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(10);

           for huy = 1:size(CoordMS,1)
                CoordQ = CoordMS(huy,:);
                CoordV4 = (sum(abs(unCorPmap(:,huy)))-1)/(size(CoordMS,1)-1);
                CoordList = CoordMS(huy+1:end,:);
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ(1,3);CoordList(hgf,3)];
                   if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                        [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                        if Rmap(hgf+huy,huy)>0
                            ySurface1 = surf(xre,yre,zre);
                            set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                            ex1_1 = 1;
                        else
                            ySurface2 = surf(xre,yre,zre);
                            set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                            ex1_2 = 1;
                        end
                   end   
                end
           end

            %%%%%%%%%%%%%%%%%%%%%%%  added to just include spheres with
            %%%%%%%%%%%%%%%%%%%%%%%  alterations
            teste = unCorPmap;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            binVet = (VetFF~=0);
            idxROIs = find(binVet);

            for i = 1:sum(VetFF~=0)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                
                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                end
                hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if exist('ex1_1') && exist('ex1_2')
                    legend([ySurface1,ySurface2],'Positive','Negative')
            else if exist('ex1_1')
                        legend(ySurface1,'Positive')
                 else
                        legend(ySurface2,'Negative')
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');

            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'WFUteste2.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'WFUteste2.nii']);
            end

            matI = aRT.dat(:,:,:);
            matI = permute(matI,[1,3,2]);
            matI = permute(matI,[2,1,3]);
            matI = permute(matI,[3,2,1]);
            matI = matI.*100;
            matI = int8(matI);
            Ds = smooth3(matI);
            Ds = int8(Ds);
            Xiso = isosurface(Ds,10);
            hiso = patch(Xiso,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3);%,'FaceVertexAlphaData',0.8,'AlphaDataMapping','Direct','CDataMapping','Direct');
            daspect([1,0.83,0.83])
            view(113,14)
            lighting gouraud
            light('Style','infinite');%'Position',[1 0 1]);
            
            saveas(GraphoFig,[Dir filesep 'Uncorrected_Graph_3D_' CovHeaders{x}],'fig')
%             saveas(GraphoFig,[Dir filesep 'Uncorrected_Graph_3D'],'png')
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_Graph_3D_' CovHeaders{x} '.png']);
            clear imgRR

            close(GraphoFig)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

%             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Color_Legend.png']);
            
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)

        end
        if get(handles.Grapho2D,'Value')  %%%%% 2D Graphos
            CoordMS2 = 3.*CoordMS;
            GraphoFig2 = figure;
            
            set(GraphoFig2,'Name','Connectome2D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                            'Color',[1 .949 .867]);
                        
            title('Uncorrected 2D Connectome','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
%             scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0 0 0]); end
            try set(gca, 'gridalpha', 1); end
            hold on
            xlim([0 3*91]);
            ylim([0 3*109]);
            zlim([0 3*91]);
            view([0 0 1])

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(15);


            for huy = 1:size(CoordMS2,1)
                CoordQ2 = CoordMS2(huy,:);
                CoordList = CoordMS2(huy+1:end,:);
                tester2 = huy<=(str_nEle+nEle)-1;
                networkNumb2 = find(tester2==1,1);
                colorV2 = cc(networkNumb2,:);
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ2(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                    if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                        if Rmap(hgf+huy,huy)>0;
                           pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                           ex2_1 = 1;
                        else
                           pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                           ex2_2 = 1;
                        end
                    end   
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%  added to just include spheres with
            %%%%%%%%%%%%%%%%%%%%%%%  alterations
            teste = unCorPmap;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            binVet = (VetFF~=0);
            idxROIs = find(binVet);

            for i = 1:sum(VetFF~=0)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS2(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.8*(5*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(5*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ(1)+(0.8*(5*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(5*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                end
                hSurface = surf(CoordV4*3.*rws+CoordQ(1,1),CoordV4*3.*res+CoordQ(1,2),CoordV4*3.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if exist('ex2_1') && exist('ex2_2')
                    legend([pLine1,pLine2],'Positive','Negative')
            else if exist('ex2_1')
                        legend(pLine1,'Positive')
                 else
                        legend(pLine2,'Negative')
                 end
            end
            
            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');

            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Contour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Contour.mat'])
            end

            daspect([1,0.95,1])
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,1);%round(size(final3D,3)/2));
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            light('Style','infinite');
            
            try
                legTmp = legend;
                IndexLegC = strfind(legTmp.String,'data1');
                IndexLeg = find(not(cellfun('isempty',IndexLegC)));
                legTmp.String(IndexLeg) = []; 
            end
            
            set(gcf, 'InvertHardCopy', 'off');
            saveas(GraphoFig2,[Dir filesep 'Uncorrected_Graph_2D_' CovHeaders{x}],'fig')
%             saveas(GraphoFig2,[Dir filesep 'Uncorrected_Graph_2D'],'png')
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_Graph_2D_' CovHeaders{x} '.png']);
            clear imgRR

            close(GraphoFig2)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gcf, 'InvertHardCopy', 'off');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])


%             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Color_Legend.png']);
            clear imgRR

            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)
        end
    else
        warndlg('No significant results with the ''uncorrected'' defined threshold ','Ops!')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDR
%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDR

if get(handles.chec2,'Value')  %%%%% FDR corrected Matrix
    FDRcor = Pmap;
    tPmap = tril(FDRcor,-1);
    tPmap = tPmap(tPmap>0);
    
    tPmap = tril(FDRcor,-1);
    tPmap = tPmap(tPmap>0);
    [~,crit_p,adj_p] = fdr_bh(tPmap,str2num(get(handles.edit2,'String')),'pdep','yes');
    stt = 1;
    res = zeros(size(g1_VAR,1),size(g1_VAR,2));
    
    for i = 1:size(g1_VAR,1)
        res(1:i,i) = zeros(i,1);
        res(i+1:end,i) = adj_p(stt:stt+size(g1_VAR,1)-i-1);
        stt = stt+(size(g1_VAR,1)-i);
    end
    resF = res+transpose(res);
    
    FDRcor = resF;

    FDRcor(FDRcor>0.05) = 1;
    FDRcor(FDRcor==0) = 1;
    FDRcor = -1.*FDRcor;
    
    if ~isequal(numel(unique(FDRcor)),1)
        fig2 = figure;
        imagesc(FDRcor)

        set(fig2,'Name','FDR-Corrected-Correlation-Matrix','Position',...
            round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                'Color',[1 .949 .867]);

        title('FDR Corrected Correlation Matrix','FontSize',14);
        xlabel('Seeds'  ,'FontSize',14);
        ylabel('Seeds','FontSize',14);
        set(gca, 'color', [0 0 0])
        try set(gca, 'gridcolor', [1 1 1]); end
        try set(gca, 'gridalpha', 1); end
        colorbar
        daspect([1 1 1])

        caxis([-1.*(str2num(get(handles.edit2,'String'))) (-1*(min(adj_p)-(0.2*min(adj_p))))])
        set(gcf, 'InvertHardCopy', 'off');
%             saveas(fig2,[Dir filesep filesep 'FDR_corrected_P_Map'],'png')
        drawnow
        imgRR = getframe(gcf);
        imwrite(imgRR.cdata, [Dir filesep 'FDR_corrected_R_Map_' CovHeaders{x} '.png']);
        clear imgRR

        saveas(fig2,[Dir filesep filesep 'FDR_corrected_R_Map_' CovHeaders{x}],'fig')
        save([Dir filesep filesep 'FDRcorRmap-VAR'],'FDRcor')
        close(fig2)

        matrixR = FDRcor;
        matrixR(matrixR==-1) = 0;
        matrixR = abs(matrixR);

            
        matrix_bin = matrixR>0;

        fidRes = fopen([DirF,filesep,'FDR-Corr_Results_Description_' CovHeaders{x} '.txt'],'w+');
        
        for r = 1:size(g1_VAR,1)
            Ridx = find(matrix_bin(r,:));
            if sum(Ridx)>0
                tester1 = r<=(str_nEle+nEle)-1;
                networkNumb1 = find(tester1==1,1);
                try
                    RnUmb1 = r - sum(nEle(1:networkNumb1-1));
                catch
                    RnUmb1 = r;
                end
                fprintf(fidRes,'ROI %d: r%d%s: \r\n',r,RnUmb1,NetLabels{networkNumb1});
                for ixR = 1:numel(Ridx)
                     if Rmap(r,Ridx(ixR))>0
                         strDir = 'Direct Correlation';
                     else
                         strDir = 'Inverse Correlation';
                     end
                    tester2 = Ridx(ixR)<=(str_nEle+nEle)-1;
                    networkNumb2 = find(tester2==1,1);
                    try
                        RnUmb2 = Ridx(ixR) - sum(nEle(1:networkNumb2-1));
                    catch
                        RnUmb2 = Ridx(ixR);
                    end
                    fprintf(fidRes,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),RnUmb2,NetLabels{networkNumb2},strDir);
                end
            end
        end
        fclose(fidRes);
        
        limiteG = (size(g1_VAR,2)/2).*1.1;

        RAIO = limiteG/2.2;
        ang=0:0.01:2*pi; 
        xp=RAIO*cos(ang);
        yp=RAIO*sin(ang);
        xpt=(limiteG/2.15)*cos(ang);
        ypt=(limiteG/2.15)*sin(ang);

        AnoRa = size(xp,2)/size(matrixR,2);
        xp2 = xp(1:AnoRa:end);
        yp2 = yp(1:AnoRa:end);

        xp2t = xpt(1:AnoRa:end);
        yp2t = ypt(1:AnoRa:end);
        xp2t = xp2t+(limiteG/2);
        yp2t = yp2t+(limiteG/2);
        xp2t = xp2t(end:-1:1);
        yp2t = yp2t(end:-1:1);

        xp2 = xp2+(limiteG/2);
        yp2 = yp2+(limiteG/2);
        yp2 = yp2(end:-1:1);
        xp2 = xp2(end:-1:1);
        
        minH = (pdist([xp2(1),yp2(1);xp2(2),yp2(2)],'euclidean'));
        
        stIndx = 1;

         for hg = 1:size(nEle,2)
             str_nEle(hg) = stIndx;
             stIndx = stIndx+nEle(hg);
         end

         cc=hsv(nOfNets);
         ccV1 = cc(1:2:end,:);
         ccV2 = cc(2:2:end,:);
         cc = [ccV1;ccV2];

        LoopAte = size(matrixR,1);
        TranspMin = 0;
        Diam = max(xp)-min(xp);
        
        circ = figure('Visible','on');
        set(circ,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                  'Color',[0 0 0]);
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Circular Connectome - DRAWING........');

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['FDR-Corr (alpha ' get(handles.edit2,'String') ')']);
        drawnow

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);
            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end
            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                    h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                    R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                    hold on
                    NewC = circ_cent(P1,P2,R);

                    Marg = (limiteG-2*RAIO)/2;
                    DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                    DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                    if DisC1>DisC2
                        centerIdx = 1;
                    else
                        centerIdx = 2;
                    end

                    AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                    AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                    if sum(abs(AngRad))>pi && any(AngRad<0)
                        idxCic = find(AngRad<0);
                        AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                    end

                    [AngRad,ksord] = sort(AngRad);

                    if isequal(ksord,[1 2])
                        TranOrd = 'descend';
                    else
                        TranOrd = 'ascend';
                    end

                    ang2 = AngRad(1):0.03:AngRad(2);
                    ang2 = [ang2,AngRad(2)];

                    xunit = R * cos(ang2) + NewC(centerIdx,1);
                    yunit = R * sin(ang2) + NewC(centerIdx,2);

                    TrnaspSca = TranspMin:(1-TranspMin)/size(xunit,2):1;
                    TrnaspSca = sort(TrnaspSca,TranOrd);
                    figure(circ)
                    for hg = 1:size(xunit,2)-1
                        h = plot(xunit(1,hg:hg+1),yunit(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                        h.Color(4) = TrnaspSca(hg);
                    end
                end
            end
        end
        set(mTextBox1,'String',['Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circ, 'InvertHardCopy', 'off');
        imgRR = getframe(circ);
        imwrite(imgRR.cdata, [DirF filesep 'FDR-Corr_Circ_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circ)

        circTS = figure('Visible','on');
        set(circTS,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                  'Color',[0 0 0],'Units','Normalized');
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        pos1 = [0.15 0.15 0.7 0.7];
        subplot('Position',pos1,'Color',[0 0 0]);
        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);
        axis off

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Circular Connectome - DRAWING........');

        TvaCrit = tinv((1-0.05)*.5,(size(g1_VAR,3)-1));

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.3  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['FDR-Corr (alpha ' get(handles.edit2,'String') ', |Tscore|>' num2str(TvaCrit) ')']);
        drawnow

        ccJet = jet(64);
        ccGre = gray(64);
        for xcv = 1:size(g1_VAR,1)
            for ycv = 1:size(g1_VAR,2)
                Tva(xcv,ycv) = tinv(1-(matrixR(xcv,ycv)*.5),(size(g1_VAR,3)-1));
            end
        end
        Tva(Tva==inf) = 0;
        TvaMax = max(max(Tva));

        cbm1 = 0;
        cbm2 = 0;
        cbm3 = 0;

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);

            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end

            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    if Rmap(i-(fact-1),vet(j))>0
                        colorVlidx = (round((Tva(i-(fact-1),vet(j))/TvaMax)*32))+32;
                        colorVl = ccJet(colorVlidx,:);
                    else
                        colorVlidx = round(((-Tva(i-(fact-1),vet(j))/TvaMax)*32))+33;
                        colorVl = ccJet(colorVlidx,:);
                    end

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                    h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                    R = ((ARCalpha^2) + (4*(h^2)))/(8*h);

                    hold on
                    NewC = circ_cent(P1,P2,R);
                    Marg = (limiteG-2*RAIO)/2;
                    DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                    DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                    if DisC1>DisC2
                        centerIdx = 1;
                    else
                        centerIdx = 2;
                    end

                    AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                    AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                    if sum(abs(AngRad))>pi && any(AngRad<0)
                        idxCic = find(AngRad<0);
                        AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                    end

                    [AngRad,ksord] = sort(AngRad);

                    if isequal(ksord,[1 2])
                        TranOrd = 'descend';
                    else
                        TranOrd = 'ascend';
                    end

                    ang2 = AngRad(1):0.03:AngRad(2);
                    ang2 = [ang2,AngRad(2)];

                    xunit = R * cos(ang2) + NewC(centerIdx,1);
                    yunit = R * sin(ang2) + NewC(centerIdx,2);
                    figure(circTS)
                    h1 = plot(xunit,yunit,'Color',colorVl,'LineWidth',3);
                    h1.Color(4) = 0.4;
                end
            end
        end

        clear A_cell Labx1
        subplot('Position',[0.04 0.01 0.9 0.015])
        imagesc(-TvaMax:(TvaMax/32):TvaMax);
        colormap(gca,'jet');
        set(gca, 'XAxisLocation', 'top','YColor',[0 0 0],'XColor',[1 1 1])
        xticks([1:64/16:64]);
        Labx1 = [(round(10*(-TvaMax:(TvaMax/8):TvaMax)))/10 TvaMax];
        A_cell = split(cellstr(num2str(Labx1)));
        xticklabels(A_cell)
        yticklabels('')
        title(gca,'T-Score Direct Corres. (hot colors) Inverse Corres. (cool colors)','Color','w')

        set(mTextBox1,'String',['Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circTS, 'InvertHardCopy', 'off');
        imgRR = getframe(circTS);
        imwrite(imgRR.cdata, [DirF filesep 'FDR-Corr_T-Score_Circ_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circTS)
        
        circ2 = figure('Visible','on');
            set(circ2,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0]);
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
        axis off
        hold on

        xlim([-1 limiteG+1]);
        ylim([-1 limiteG+1]);

        daspect([1 1 1])
        cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

         mTextBox1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.3  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox1,'String','Organic Circular Connectome - DRAWING........');

         mTextBox2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
             'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
             'ForegroundColor',[1 .949 .867],'FontWeight','bold');
         set(mTextBox2,'String',['FDR-Corr (alpha ' get(handles.edit2,'String') ')']);
        drawnow

        daspect([1 1 1])
        RAIO2 = RAIO*0.15;        
        xp3=RAIO2*cos(ang);
        yp3=RAIO2*sin(ang);

        xp3 = xp3(1:AnoRa:end);
        yp3 = yp3(1:AnoRa:end);
        xp3 = flip(xp3+(limiteG/2));
        yp3 = flip(yp3+(limiteG/2));

        RAIO3 = RAIO*0.85;        
        xp4=RAIO3*cos(ang);
        yp4=RAIO3*sin(ang);

        idxNR = ceil((str_nEle+(nEle./2))-1);
        xp4 = flip(xp4(1:AnoRa:end));
        yp4 = flip(yp4(1:AnoRa:end));
        xp4 = xp4(idxNR);
        yp4 = yp4(idxNR);

        xp4 = xp4+(limiteG/2);
        yp4 = yp4+(limiteG/2);

        RAIO4 = RAIO*0.95;        
        xp5=RAIO4*cos(ang);
        yp5=RAIO4*sin(ang);

        xp5 = xp5((AnoRa/2):AnoRa:end);
        yp5 = yp5((AnoRa/2):AnoRa:end);
        xp5 = flip(xp5+(limiteG/2));
        yp5 = flip(yp5+(limiteG/2));

        for i = fact:LoopAte
            drawnow
            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);
            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end
            vet = matrixR(i-(fact-1),:);
            vet = find(vet);

            if ~isempty(vet)
                bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                if abs(Angle)>pi/2
                    Angle = Angle - pi;
                    set(bt,'HorizontalAlignment','right')
                end
                set(bt,'Rotation',rad2deg(Angle))
                for j = 1:size(vet,2)

                    P1 = [xp2(i),yp2(i)];
                    P2 = [xp2(vet(j)),yp2(vet(j))];

                    testerj = vet(j)<=(str_nEle+nEle)-1;
                    networkNumbj = find(testerj==1,1);
                    colorVj = cc(networkNumbj,:);
                    xxF = [];
                    yyF = [];
                    if networkNumb~= networkNumbj
                        opt1 =  round(i-(i-vet(j))/2);
                        opt2 = round(opt1+(LoopAte/2));
                        if opt2<=LoopAte
                            dopt1 = pdist([xp3(opt1),yp3(opt1);P1(1),P1(2)]);
                            dopt2 = pdist([xp3(opt2),yp3(opt2);P1(1),P1(2)]);
                            if dopt1>dopt2
                                dopt = opt2;
                            else
                                dopt = opt1;
                            end
                        else
                            dopt = opt1;
                        end

                        [xx,yy]=fillline([P1(1),P1(2)],[xp4(networkNumb),yp4(networkNumb)],15);
                        [xx2,yy2]=fillline([xp4(networkNumb),yp4(networkNumb)],[xp3(dopt),yp3(dopt)],20);
                        [xx3,yy3]=fillline([xp3(dopt),yp3(dopt)],[xp4(networkNumbj),yp4(networkNumbj)],20);
                        [xx4,yy4]=fillline([xp4(networkNumbj),yp4(networkNumbj)],[P2(1),P2(2)],15);

                        xxF = [xx,xx2,xx3,xx4];
                        yyF = [yy,yy2,yy3,yy4];
                        xxF = [xxF(1);smoothLNI(xxF(2:end-1),30);xxF(end)]';
                        yyF = [yyF(1);smoothLNI(yyF(2:end-1),30);yyF(end)]';
                        TranOrd = 'descend';

                        TrnaspSca = [zeros(1,size(xxF,2)*0.3)+TranspMin,TranspMin:(1-TranspMin)/(size(xxF,2)*0.7):1];
                        TrnaspSca = sort(TrnaspSca,TranOrd);
                        figure(circ2)
                        if networkNumb~=networkNumbj
                            for hg = 1:(size(xxF,2)-1)*0.7
                                h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                h.Color(4) = TrnaspSca(hg);
                            end
                        else
                            if i<vet(j)
                                h = plot(xxF,yyF,'Color',colorV,'LineWidth',1);
                            end
                        end

                    else

                        ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                        h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                        R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                        NewC = circ_cent(P1,P2,R);
                        Marg = (limiteG-2*RAIO)/2;
                        DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                        DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                        if DisC1>DisC2
                            centerIdx = 1;
                        else
                            centerIdx = 2;
                        end

                        AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                        AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                        if sum(abs(AngRad))>pi && any(AngRad<0)
                            idxCic = find(AngRad<0);
                            AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                        end

                        ang2 = AngRad(1):0.03:AngRad(2);
                        ang2 = [ang2,AngRad(2)];

                        xunit = R * cos(ang2) + NewC(centerIdx,1);
                        yunit = R * sin(ang2) + NewC(centerIdx,2);

                        h = plot(xunit,yunit,'Color',colorV.*0.5,'LineWidth',2);
                    end

                    if vet(j)<10
                        strf = ['0' num2str(vet(j))];
                    else
                        strf = num2str(vet(j));
                    end
                    bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                        colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                        'VerticalAlignment','middle');

                    Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                end
            end
        end

        set(mTextBox1,'String',['Organic Circular Connectome (' CovHeaders{x} ')']);
        drawnow

        set(circ2, 'InvertHardCopy', 'off');
        imgRR = getframe(circ2);
        imwrite(imgRR.cdata, [DirF filesep 'FDR-Corr_OrganicCirc_Map_' CovHeaders{x} '.png']);
        clear imgRR
        close(circ2)
        
        if get(handles.Graphos3D,'Value') %%%%% FDR corrected 3D Grapho
            GraphoFig3 = figure;
            set(GraphoFig3,'Name','Connectome3D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                            'Color',[1 .949 .867]);
                        
            title('FDR Corrected 3D Connectome','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0.8 0.8 0.8]); end
            try set(gca, 'gridalpha', 1); end
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            hold on
            view([-1 -1 1])

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(15);
             for huy = 1:size(CoordMS,1)
                CoordQ = CoordMS(huy,:);
                CoordList = CoordMS(huy+1:end,:);
                CoordV4 = sum(FDRcor(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ(1,3);CoordList(hgf,3)];
                    if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                        [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
                        if Rmap(hgf+huy,huy)>0;
                            ySurface1 = surf(xre,yre,zre);
                            set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                            ex3_1 = 1;
                        else
                            ySurface2 = surf(xre,yre,zre);
                            set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                            ex3_2 = 1;
                        end
                    end   
                end
             end
            
            teste = FDRcor;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            binVet = (VetFF~=0);
            idxROIs = find(binVet);

            for i = 1:sum(VetFF~=0)

                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                end
                hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end

            if exist('ex3_1') && exist('ex3_2')
                    legend([ySurface1,ySurface2],'Positive','Negative')
            else if exist('ex3_1')
                        legend(ySurface1,'Positive')
                 else
                        legend(ySurface2,'Negative')
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'WFUteste2.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'WFUteste2.nii']);
            end

            matI = aRT.dat(:,:,:);
            matI = permute(matI,[1,3,2]);
            matI = permute(matI,[2,1,3]);
            matI = permute(matI,[3,2,1]);
            matI = matI.*100;
            matI = int8(matI);
            Ds = smooth3(matI);
            Ds = int8(Ds);
            Xiso = isosurface(Ds,10);
            hiso = patch(Xiso,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3);%,'FaceVertexAlphaData',0.8,'AlphaDataMapping','Direct','CDataMapping','Direct');
            daspect([1,0.83,0.83])
            view(113,14)
            lighting gouraud

            light('Style','infinite');%'Position',[1 0 1]);
            
            saveas(GraphoFig3,[Dir filesep 'FDR_Corrected_Graph_3D_' CovHeaders{x}],'fig')
            set(gcf, 'InvertHardCopy', 'off');
            
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'FDR_Corrected_Graph_3D_' CovHeaders{x},'.png']);
            clear imgRR

            close(GraphoFig3)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Color_Legend.png']);
            clear imgRR

            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)
            
        end
        
        if get(handles.Grapho2D,'Value')  %%%%% FDR corrected 2D Grapho
            CoordMS2 = 3.*CoordMS;
            GraphoFig4 = figure;
            set(GraphoFig4,'Name','Connectome2D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                            'Color',[1 .949 .867]);
            title('FDR Corrected 2D Connectome','FontSize',14);
            
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0 0 0]); end
            try set(gca, 'gridalpha', 1); end
            hold on
            xlim([0 3*91]);
            ylim([0 3*109]);
            zlim([0 3*91]);
            view([0 0 1])
            
             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

            cc=hsv(nOfNets);
            [rws,res,rdw] = sphere(15);
            
            for huy = 1:size(CoordMS2,1)
                CoordQ2 = CoordMS2(huy,:);
                CoordList = CoordMS2(huy+1:end,:);
                CoordV4 = sum(FDRcor(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS2,1));
                tester2 = huy<=(str_nEle+nEle)-1;
                networkNumb2 = find(tester2==1,1);
                colorV2 = cc(networkNumb2,:);
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ2(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                    if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                        if Rmap(hgf+huy,huy) > 0
                           pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                           ex4_1 = 1;
                        else
                           pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                           ex4_2 = 1;
                        end
                    end   
                end
            end
            
            teste = FDRcor;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            binVet = (VetFF~=0);
            idxROIs = find(binVet);

            for i = 1:sum(VetFF~=0)

                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS2(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);
                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.9*(5*CoordV4)),CoordQ(2)+(0.9*(5*CoordV4)),CoordQ(3)+(0.9*(5*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',11);
                else
                    text(CoordQ(1)+(0.9*(5*CoordV4)),CoordQ(2)+(0.9*(5*CoordV4)),CoordQ(3)+(0.9*(5*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',11);
                end
                hSurface = surf(CoordV4.*rws.*3+CoordQ(1,1),CoordV4.*3.*res+CoordQ(1,2),CoordV4.*3.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            if exist('ex4_1') && exist('ex4_2')
                    legend([pLine1,pLine2],'Positive','Negative')
            else if exist('ex4_1')
                        legend(pLine1,'Positive')
                 else
                        legend(pLine2,'Negative')
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');

            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Contour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Contour.mat'])
            end

            daspect([1,0.95,1])
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,1);%round(size(final3D,3)/2));
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            light('Style','infinite');
            
            set(gcf, 'InvertHardCopy', 'off');
            saveas(GraphoFig4,[Dir filesep 'FDR_Corrected_Graph_2D_' CovHeaders{x}],'fig')
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'FDR_Corrected_Graph_2D_' CovHeaders{x},'.png']);
            clear imgRR

            close(GraphoFig4)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            drawnow
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [Dir filesep 'Color_Legend.png']);
            clear imgRR

            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)
            
        end
    else
        warndlg('No significant results with FDR correction','Ops!')
    end
end
clear('ex1_1','ex2_1','ex1_2','ex2_2','ex3_1','ex3_2','ex4_2','ex4_1')
end

set(handles.text13,'String','Done!')
disp('------------')
disp('You can open the *.fig files on the Matlab to interact (rotate, zoom, edit...) with the figures!')

function ROIfiles_Callback(hObject, eventdata, handles)
if get(handles.ROIfiles,'Value')
    set(handles.CoodList,'Value',0)
    set(handles.addNODES,'String','Add ROIs images')
else
    set(handles.CoodList,'Value',1)
    set(handles.addNODES,'String','Add Coordinate list')
end

function CoodList_Callback(hObject, eventdata, handles)
if get(handles.CoodList,'Value')
    set(handles.ROIfiles,'Value',0)
    set(handles.addNODES,'String','Add Coordinate list')
else
    set(handles.ROIfiles,'Value',1)
    set(handles.addNODES,'String','Add ROIs images')
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
