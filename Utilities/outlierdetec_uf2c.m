function [OutInd,Cutoff] = outlierdetec_uf2c(input,OutType,direction)
%
% OutInd = outlierdetec_uf2c(input,OutType,direction)
%
% INPUTS
% input: Numerical vector
%
% OutType: the Outlier type
%          'major' (only very outliers: 3*(Q3-Q1))
%          'minor' (more easy to be an outliers: 1.5*(Q3-Q1)) **STANDARD OPTION IF EMPTY**
%          'severe' (very easy to be considered an outliers: 1*(Q3-Q1))
%
% direction: define if you want only upper outlies, lower outliers or both direction
%           'twosides' (both direction outliers) **STANDARD OPTION IF EMPTY**
%           'upperside' (only outliers with hight values)
%           'lowerside' (only outliers with low values)
%
%  OUTPUT: The indices of the outliers on the input vector

% UF²C - User Friendly Functional Connectivity
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
    
    if ~exist('direction','var')
        direction = 'twosides';
    end
    if ~exist('OutType','var')
        OutType = 'minor';
    end

    if ~isequal(OutType,'major') && ~isequal(OutType,'minor') && ~isequal(OutType,'severe')
        error('The Outlier type should be ''major'', ''minor'' or ''severe''')
        return
    end

    if ~isequal(direction,'twosides') && ~isequal(direction,'upperside') && ~isequal(direction,'lowerside')
        error('The Outlier Direction should be ''twosides'', ''upperside'' or ''lowerside''')
        return
    end
    
    vet2 = sort(input);
    sizeQ = numel(vet2)/4;

    if isequal(rem(sizeQ,1),0)
        Q1 = (vet2(floor(sizeQ))+vet2(floor(sizeQ)+1))/2;
        Q3 = (vet2(floor(3*sizeQ))+vet2(floor(3*sizeQ)+1))/2;
    else
        Q1 = vet2(ceil(sizeQ));
        Q3 = vet2(ceil(3*sizeQ));
    end

    IntQRng = Q3-Q1;

    switch OutType
        case 'major'
            InnerFacmaj = IntQRng*3;
            InnerQ1maj = Q1-InnerFacmaj;
            InnerQ3maj = Q3+InnerFacmaj;
            switch direction
                case 'twosides'
                    OutInd = find(input<InnerQ1maj | input>InnerQ3maj);
                    Cutoff = [InnerQ3maj;InnerQ1maj];
                case 'upperside'
                    OutInd = find(input>InnerQ3maj);
                    Cutoff = InnerQ3maj;
                case 'lowerside'
                    OutInd = find(input<InnerQ1maj);
                    Cutoff = InnerQ1maj;
            end
        case 'minor'
            InnerFacmin = IntQRng*1.5;
            InnerQ1min = Q1-InnerFacmin;
            InnerQ3min = Q3+InnerFacmin;
            switch direction
                case 'twosides'
                    OutInd = find(input<InnerQ1min | input>InnerQ3min);
                    Cutoff = [InnerQ3min;InnerQ1min];
                case 'upperside'
                    OutInd = find(input>InnerQ3min);
                    Cutoff = InnerQ3min;
                case 'lowerside'
                    OutInd = find(input<InnerQ1min);
                    Cutoff = InnerQ1min;
            end
        case 'severe'
            InnerFacsevere = IntQRng*1;
            InnerQ1severe = Q1-InnerFacsevere;
            InnerQ3severe = Q3+InnerFacsevere;
            switch direction
                case 'twosides'
                    OutInd = find(input<InnerQ1severe | input>InnerQ3severe);
                    Cutoff = [InnerQ3severe;InnerQ1severe];
                case 'upperside'
                    OutInd = find(input>InnerQ3severe);
                    Cutoff = InnerQ3severe;
                case 'lowerside'
                    OutInd = find(input<InnerQ1severe);
                    Cutoff = InnerQ1severe;
            end
    end
end


