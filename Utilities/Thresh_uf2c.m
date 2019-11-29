function ThreInpu = Thresh_uf2c(image,type,QVthre,CluSiz)
%
% Brunno Machado de Campos
% University of Campinas, 2017
%
% Copyright (c) 2017, Brunno Machado de Campos
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
%
% INPUTS: 
%   image: can be:
%        *.mat file with a matrix variable (2D, 3D or 4D)
%        NIfTI (*.nii) file of an image (3D or 4D)
%        Existent Workspace matrix variable (2D, 3D or 4D)
%
%   type: The input 'image' type. The possible values are:
%        'mfile' (for *.mat file with variables)
%        'nii' (for NIfTI (*.nii) image)
%        'var' (existent Workspace variable)
%
%   QVthre: A absolute threshold value. The possible values are:
%        Double precision value. E.g. 0.2 or 0 or 100 (depending of the matrix values and your interest)
%        Empty value: []. In this case the image will be not thresholded
%
%   CluSiz: The minimal cluster size.
%         Double precision value. E.g. 50 or 10 or 5....
%         Value of zero: 0. In this case the image will be not thresholded by cluster size
% 
% 
% OUTPUT: The outputs will vary acoording to the input type:
%   'mfile': The output will be a string with the whole filename of 
%         the new modified mat-file saved in the input mat-file folder. 
%         The variable name inside the mat-file will not be changed.
%
%   'nii': The output will be a string with the whole filename of 
%         the new modified NIfTI file saved in the input image folder.
%         
%   'var': A new variable with the defined name and with the 
%         transformed matrix. E.g. NewVar = Thresh_uf2c(Map,'var',0.2,50);
% 
% Examples: 
%   ThreInpu = Thresh_uf2c('C:\Users\spmT_0001.nii','nii',3.2,10);
%
%   ThreInpu = Thresh_uf2c('C:\Users\3DMap.mat','mfile',[],50);
%   ThreInpu = Thresh_uf2c('C:\Users\3DMap.mat','mfile',0.3,0);
%
%   ThreInpu = Thresh_uf2c(Map,'var',0.2,50);
% 

    if ~strcmp(type,'mfile') && ~strcmp(type,'nii') && ~strcmp(type,'var')
        error('ErrorTests:convertTest',...
            ['Error using ExtThre_uf2c. No valid input image ''type'' specified.\n',...
            'The options are: \n',...
                '\t''mfile'' (*.mat file with variables)\n',...
                '\t''nii'' (NIfTI image, *.nii)\n',...
                '\t''var'' (Existent Workspace variable)\n'])
        return
    end

    switch type
        case 'mfile'
            a = load(image);
            [path,file] = fileparts(image);
            value = fields(a);
            value = value{1};
            eval(sprintf('%s = a.%s;',value,value));
            tester = eval(value);
            if isequal(numel(size(tester)),3) || isequal(numel(size(tester)),2)
                intLoop = 1;
            else
                intLoop = size(tester,4);
            end
        case 'nii'
            [path,file] = fileparts(image);
            stru = nifti(image);
            if isequal(numel(stru.dat.dim),3)
                tester = stru.dat(:,:,:);
                intLoop = 1;
            end
            if isequal(numel(stru.dat.dim),4)
                tester = stru.dat(:,:,:,:);
                intLoop = stru.dat.dim(4);
            end
        case 'var'
            tester = image;
            if isequal(numel(size(tester)),3) || isequal(numel(size(tester)),2)
                intLoop = 1;
            else
                intLoop = size(tester,4);
            end
    end
    
    if ~isempty(QVthre)
        tester2 = tester>QVthre;
    else
        tester2 = tester;
    end

    if isequal(intLoop,1)
        L = bwlabeln(tester2);
        vet = unique(L);
        vet = vet(2:end);
        Fmat = zeros(size(tester2,1),size(tester2,2),size(tester2,3));
        nClu = 0;
        for i = 1:size(vet,1)
            Mat = double(L==i);
            if sum(sum(sum(Mat)))>CluSiz
               Fmat = Fmat + Mat;
               nClu = nClu +1;
            end
        end
        FmatF = Fmat.*tester;
    else
        FmatF = zeros(size(tester2));
        for gf = 1:intLoop
            L = bwlabeln(tester2(:,:,:,gf));
            vet = unique(L);
            vet = vet(2:end);
            Fmat = zeros(size(tester2,1),size(tester2,2),size(tester2,3));
            nClu = 0;
            for i = 1:size(vet,1)
                Mat = double(L==i);
                if sum(sum(sum(Mat)))>CluSiz
                   Fmat = Fmat + Mat;
                   nClu = nClu +1;
                end
            end
            FmatF(:,:,:,gf) = Fmat.*tester(:,:,:,gf);
            clear L vet
        end
    end

    switch type
        case 'mfile'
            eval(sprintf('%s = FmatF;',value))
            save([path filesep 'Vthre' num2str(QVthre) '_Kmin' num2str(CluSiz) '_' file '.mat'],sprintf('%s',value),'-v7.3');            ThreInpu = [path filesep 'Vthre' num2str(QVthre) '_Kmin' num2str(CluSiz) '_' file];
            ThreInpu = [path filesep 'Vthre' num2str(QVthre) '_Kmin' num2str(CluSiz) '_' file '.mat'];
        case 'nii'
            stru2 = stru;
            stru2.dat.fname = [path filesep 'Vthre' num2str(QVthre) '_Kmin' num2str(CluSiz) '_' file '.nii'];
            if isequal(intLoop,1)
                stru2.dat(:,:,:) = FmatF;
            else
                stru2.dat(:,:,:,:) = FmatF;
            end
            create(stru2)
            ThreInpu = [path filesep 'Vthre' num2str(QVthre) '_Kmin' num2str(CluSiz) '_' file '.nii'];
        case 'var'
            ThreInpu = FmatF;
    end
end   
    
    
    
    
