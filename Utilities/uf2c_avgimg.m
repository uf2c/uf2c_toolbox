function [AvgImage,AvgImageMASK] = uf2c_avgimg(images,thresh,ttype,OutDir,OutName,SaveNiis)
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
%
%
% INPUTS: 
%   images [OPTIONAL]: Cell array of strings with the images filenames.
%       E.g.: {'E:\Images\CASE01.nii';'E:\Images\CASE02.nii';'E:\Images\CASE03.nii'}
% 
%   thresh [OPTIONAL]: Double variable. Absolute threshold to convert voxels into zero. 
%       E.g.: 0.2
% 
%   ttype [OPTIONAL]: String variable. 
%   OPTIONS: 
%   'PRE' The threshold will be applied to each image before the average.
%   'POS' The threshold will be applied only to the average image.
%   'PREPOS' The threshold will be applied to each image before the
%   average and to the final average image;
% 
%   OutDir [OPTIONAL]: String variable with the output path
%       E.g.: 'E:\Images\'
%
%   OutNam [OPTIONAL]: String variable with the output average image filename
%       E.g.: 'MyAverage.nii'
%   
%   SaveNiis [OPTIONAL]: Logical variable
%   OPTIONS: 
%   0 No average image NIfTI file will be saved
%   1 The average image NIfTI file will be saved
%
% OUTPUTS:
%   The average image in NIfTI format
%
%   Binary average image in NIfTI format
%
%   AvgImage: The average image matrix
% 
%   AvgImageMASK: The binary average image matrix

if ~exist('images','VAR')
    [images,pathG1] = uigetfile({'*.nii','NIfTI files'},...
        'Add all images to compute the average','MultiSelect','on');
else
    if ~iscell(images)
        warndlg('You should add a cell array of strings with images filenames','Attention.')
        return
    end
end
if ~exist('images','VAR')
    thresh = 'none';
end
if ~exist('thresh','VAR') || isempty(thresh)
    thresh = 'none';
end

if ~exist('SaveNiis','VAR')
    SaveNiis = 1;
end

if SaveNiis
    if ~exist('OutDir','VAR')
        if exist('pathG1','VAR')
            OutDir = pathG1;
        else
            OutDir = pwd;
        end
    end
    if ~exist('OutName','VAR')
        OutName = 'AverageImage.nii';
    end
end

vol4D = 1;

for i = 1:numel(images)
    Strutmp = nifti(images{i});
    MatF(:,:,:,i) = Strutmp.dat(:,:,:,vol4D);
end

if isequal(thresh,'none')
    AvgImage = sum(MatF,4)./size(MatF,4);
else
    if ~exist('ttype','VAR')
        ttype = 'PRE';
    end
    switch ttype
        case 'PRE'
            MatF(MatF < thresh) = 0;
            AvgImage = sum(MatF,4)./size(MatF,4);
        case 'POS'
            AvgImage = sum(MatF,4)./size(MatF,4);
            AvgImage(AvgImage < thresh) = 0;
        case 'PREPOS'
            MatF(MatF < thresh) = 0;
            AvgImage = sum(MatF,4)./size(MatF,4);
            AvgImage(AvgImage < thresh) = 0;
    end
end

AvgImageMASK = AvgImage;
AvgImageMASK(AvgImageMASK > 0) = 1;

if SaveNiis
    StruF = Strutmp;
    StruF.dat.fname = [OutDir filesep OutName];
    StruF.dat.dim = size(AvgImage);
    StruF.dat(:,:,:) = AvgImage;
    create(StruF)

    StruFbin = Strutmp;
    StruFbin.dat.fname = [OutDir filesep OutName(1:end-4) 'BinMask.nii'];
    StruFbin.dat.dim = size(AvgImageMASK);
    StruFbin.dat(:,:,:) = AvgImageMASK;
    create(StruFbin)
end
