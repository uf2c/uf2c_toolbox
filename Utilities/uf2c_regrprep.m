function [Frgval,Regstr,figu,thicksL] = uf2c_regrprep(nOfSub,checkMOV,pathFunc,dirname,fileFunc,checkOsc,checkReg,ScreSize,MeanWM,MeanCSF,dyn4,nr)
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2018
%
% Copyright (c) 2018, Brunno Machado de Campos
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

% This function mask the functional data and extract globals signals for
% regressions.
% the function assumes that the "dirname" folder contains normalized
% segmented mpas of the structural image (fileStru): wc1*fileStru, 
% wc2*fileStru, wc3*fileStru.
%
% The imputs are:
%
% checkMOV: binary variable indicating if the regression for movement will be performed (UF²C default: 1)
% pathFunc: directory where the functional image folder is/ string
% dirname: name of the folder where the functional image is/ string
% fileFunc: name of the functional image (no path)/ string
% checkOsc: binary variable indicating if the global signals regression will be performed (UF²C default: 1)
% checkReg: binary variable indicating if any regression will be performed (UF²C default: 1)
% ScreSize: vector (1x2) with the screen resolution
% MeanWM: vector (3x4th dim) with the PCA (3 firsts) white matter signals time
%   series
% MeanCSF: vector (3x4th dim) with the PCA (3 firsts) CSF signals time
%   series
% dyn4: 4th dim size (double variable)
%
% Ex: [Frgval,Regstr] = regrprep_uf2c(1,'C:\Users\fMRI\Tests','Vol1','Vol1_fMRI.nii',1,1,[1600,1200],MeanWM,MeanCSF,180)
%
% The outputs are:
%
% Frgval: matrix containning the regressions
% Regstr: string describing what the matrix contains
%

    r1 = '';
    r2 = '';
    r3 = '';

    if isequal(checkMOV,1)
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Preparing regressors\r\n');

        try
            rpfile = textread([pathFunc,dirname,filesep,'rp_',fileFunc(1:end-3),'txt']); % read alignment file output
        catch
            rpfile = textread(fileFunc); % read alignment file output
        end
        rgval1 = rpfile; %add movement parameters to regress matrix

        Min_rp_1 = min(rgval1(1:end,1));
        Min_rp_2 = min(rgval1(1:end,2));
        Min_rp_3 = min(rgval1(1:end,3));
        Min_rp_4 = min(rgval1(1:end,4));
        Min_rp_5 = min(rgval1(1:end,5));
        Min_rp_6 = min(rgval1(1:end,6));

        rgval1(1,:) = 0;

        rgval1(1:end,1) = rgval1(1:end,1) + abs(Min_rp_1);
        rgval1(1:end,2) = rgval1(1:end,2) + abs(Min_rp_2);
        rgval1(1:end,3) = rgval1(1:end,3) + abs(Min_rp_3);
        rgval1(1:end,4) = rgval1(1:end,4) + abs(Min_rp_4);
        rgval1(1:end,5) = rgval1(1:end,5) + abs(Min_rp_5);
        rgval1(1:end,6) = rgval1(1:end,6) + abs(Min_rp_6);

        r1 = 'mov-';
    end

    clear Min_rp_1 Min_rp_2 Min_rp_3 Min_rp_4 Min_rp_5 Min_rp_6

    if isequal(checkOsc,1)
%         rgval2 = transpose([MeanWM;MeanCSF]); %add mean white matter and CSF parameter to regress matrix
        rgval2 = [MeanWM,MeanCSF]; %add PCAs of white matter and CSF parameter to regress matrix
        rgval2Min = min(rgval2);
        
        for i = 1:size(rgval2,2)
            rgval2(:,i) = rgval2(:,i) + abs(rgval2Min(i));
        end
        rgval2Max = max(rgval2);
        for i = 1:size(rgval2,2)
            rgval2(:,i) = rgval2(:,i)./abs(rgval2Max(i));
        end
        for i = 1:size(MeanCSF,2)
            if ~isempty(MeanWM) % condition for GeneralPreproc regressors
                OSCthicksLWM{i,1} = ['WM PC',num2str(i)];
            else
                OSCthicksLWM = [];
            end
            OSCthicksLCSF{i,1} = ['CSF PC',num2str(i)];
        end
        OSCthicksL = [OSCthicksLWM;OSCthicksLCSF];
        r2 = 'osc-';
    end
    
    UF2Cdir = which('uf2c');
    tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];
    numofVoluF = size(MeanWM,1);
    if isequal(checkReg,1)
        if ~isequal(nr,0)
            rgval3 = load([tmpDIR 'additional_Reg.mat']);
            nCol = 1;
            for yt = 1:numofVoluF:nr*numofVoluF
                rgval3r(:,nCol) = rgval3.regVets(yt:yt+(numofVoluF-1),nOfSub);
                nCol = nCol +1;
            end
            rgval3rMin = min(rgval3r);

            for i = 1:size(rgval3r,2)
                rgval3r(:,i) = rgval3r(:,i) + abs(rgval3rMin(i));
            end
            rgval3rMax = max(rgval3r);
            for i = 1:size(rgval3r,2)
                rgval3r(:,i) = rgval3r(:,i)./abs(rgval3rMax(i));
            end
            r3 = 'addt';
            for hgf = 1:nCol-1
                addRegTick{hgf,1} = ['Addit.Reg' num2str(hgf)];
            end
        end
    end

    if isequal(r3,'addt')
        if size(unique(rgval3r),1)==1
            r3 = '';
        end
    end

    Regstr = [r1,r2,r3];

    if isequal(Regstr,'osc-')
        rgval = rgval2;
        thicksL = [OSCthicksL;{'Ones'}];
    end
    if isequal(Regstr,'mov-')
        rgval = rgval1;
        thicksL = {'Displ. X';'Displ. Y';'Displ. Z';'Rot. Pitch';'Rot. Row';'Rot. Yaw';'Ones'};
    end
    if isequal(Regstr,'addt')
        rgval = rgval3r;
        thicksL = {'Addit. reg';'Ones'};
    end
    if isequal(Regstr,'mov-osc-')
        rgval = [rgval2,rgval1];
        thicksL = [OSCthicksL;{'Displ. X';'Displ. Y';'Displ. Z';'Rot. Pitch';'Rot. Row';'Rot. Yaw';'Ones'}];
    end
    if isequal(Regstr,'mov-addt')
        rgval = [rgval1,rgval3r];
        thicksL = [{'Displ. X';'Displ. Y';'Displ. Z';'Rot. Pitch';'Rot. Row';'Rot. Yaw'};addRegTick;{'Ones'}];
    end
    if isequal(Regstr,'mov-osc-addt')
        rgval = [rgval3r,rgval2,rgval1];
        thicksL = [addRegTick;OSCthicksL;{'Displ. X';'Displ. Y';'Displ. Z';'Rot. Pitch';'Rot. Row';'Rot. Yaw'};{'Ones'}];
    end
    if isequal(Regstr,'osc-addt')
        rgval = [rgval3r,rgval2];
        thicksL = [addRegTick;OSCthicksL;{'Ones'}];
    end
    if isequal(Regstr,'')
       fprintf('The use of regressos are strongly recommended! Continuing...\n');
       Frgval = [];
       Regstr = [];
       figu = [];
    end

    if ~isequal(Regstr,'')
        for rte = 1:size(rgval,2) % all regressors are normalized: minimal value = 0 maximal value = 1. the first line is a copy of the second.
            MinV = min(rgval(2:end,rte));
            MaxV = max(rgval(2:end,rte));
            for ttt = 1:dyn4
                Frgval(ttt,rte) = rgval(ttt,rte)./MaxV;
            end
            Frgval(:,rte) = Frgval(:,rte) - mean(Frgval(:,rte)); % Mean Centering
        end

        Frgval = [Frgval,ones(dyn4,1)];
        fprintf('Done! ============== %s\r\n',spm('time'));    
    end
    
    figu = figure('Visible','on');
    set(figu,'Name','RegressionMatrix',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.5 ScreSize(1)*.3]),...
        'Color',[1 0.94 0.86]);
    imagesc(Frgval);
    grid('off')

    colormap(gray)
    title('Regression Matrix','FontSize', 14);
    set(gca,'XTick',1:size(thicksL,1))
    set(gca,'XTickLabel',thicksL)
    xticklabel_rotate([],90,[],'Fontsize',10)
    drawnow
end