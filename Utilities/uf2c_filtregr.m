function finalEPI = uf2c_filtregr(finalEPI,BinFazFilt,Dx,Dy,Dz,Dt,Frgval,HdLOW,HdHIGH,aFil,tmp1,tmp2,Regstr)

    finalEPIresh = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';

    testerMatx = finalEPIresh;
    testerMatx(testerMatx~=0) = 1;
    vetBin = sum(testerMatx,1);
    vetBin(vetBin<Dt) = 0; % Removes from the mask, voxels with time series with more than 5% of zero values
    vetBin(vetBin==Dt) = 1;
    vetBinMat = bsxfun(@and,vetBin,ones(1,Dt)');
    finalEPIresh = finalEPIresh.*vetBinMat;

    finalEPIresh(:,~any(finalEPIresh,1) ) = [];  %columns
    finalEPIreshTMP = zeros(size(finalEPIresh));
    
    clear testerMatx vetBinMat

%     Dyns = 180;
%     TR = 2;
%     FreqLowP = 0.1;
%     FreqHighP = 0.008;
%     FreqStep = 1/Dyns;
%     PosSteps = 0:FreqStep:0.5-FreqStep;
%     [~,DiffLP] = min(abs(PosSteps-FreqLowP));
%     RealFreqLowP = PosSteps(DiffLP);
%     [~,DiffHP] = min(abs(PosSteps-FreqHighP));
%     RealFreqHighP = PosSteps(DiffHP);
%     RemFreqs = PosSteps(DiffHP:DiffLP);
%     fs = 1/TR;
%     t = (0:TR:Dyns*TR)';
%     for i = 1:numel(RemFreqs)
%         data(:,i) = sin(2*pi*RemFreqs(i)*t)';
%     end
%     data = data(1:end-1,:);
% 
%     Frgval = [Frgval,data];
    
    if ~isempty(gcp('nocreate'))
        multiWaitbar('Filtering & Regression', 0.5);
    end
    
    NumT1 = Dt; % number of time points
    if ~isequal(Regstr,'') && isequal(BinFazFilt,1)
        if ~isempty(HdLOW)
            HdLOW_Numerator = HdLOW.Numerator;
            fLPF = 1;
        else
            fLPF = 0;
        end
        if ~isempty(HdHIGH)
            HdHIGH_Numerator = HdHIGH.Numerator;
            fHPF = 1;
        else
            fHPF = 0;
        end
        parfor tt3 = 1:size(finalEPIresh,2)

            multiWaitbar('Filtering & Regression', tt3/size(finalEPIresh,2));
            SignCorr = finalEPIresh(:,tt3);
            MMM = mean(SignCorr);

            % Detrend
            SignCorr = detrend(SignCorr,'linear');  %detrend with breakpoints over 128 seconds
            SignCorr = SignCorr + MMM;
            MMM = mean(SignCorr);

            % Regression. basically includes: Movement regressors(6 parameters), mean White Matter signal and mean CSF signal
            b = regress(SignCorr,Frgval); 
            WReg  = Frgval*b;
            s_corr = SignCorr - WReg;
            s_corr = s_corr + MMM;
            
            % Low Pass Filtering
            if fLPF % Case of no low pass filter
                if NumT1<=tmp1 % if the n of dynamics if lower than 3 times the size of the filtering object
                    s_filt2 = repmat(s_corr,ceil(tmp1/NumT1),1); % we use the center copy, avoinding phase distortion
                else
                    s_filt2 = repmat(s_corr,3,1);
                end
                s_filtAB1LOW = filtfilt(HdLOW_Numerator,aFil,s_filt2);
                s_filtAB1LOW = s_filtAB1LOW((1+NumT1):(2*NumT1));
            else
                s_filtAB1LOW = s_corr;
            end

            % High Pass Filtering
             if fHPF % Case of no high pass filter
                if NumT1<=tmp2
                    s_filtAB1LOW = repmat(s_filtAB1LOW,ceil(tmp2/NumT1),1);
                else
                    s_filtAB1LOW = repmat(s_filtAB1LOW,3,1);
                end
                s_filtAB1LOW_HIGH = filtfilt(HdHIGH_Numerator,aFil,s_filtAB1LOW);
                s_filtAB1LOW_HIGH = s_filtAB1LOW_HIGH((1+NumT1):(2*NumT1));
             else
                 s_filtAB1LOW_HIGH = s_filtAB1LOW;
             end

             finalEPIreshTMP(:,tt3) = s_filtAB1LOW_HIGH(:) + b(end);

        end
        finalEPIresh = finalEPIreshTMP;
    end
    if ~isequal(Regstr,'') && isequal(BinFazFilt,0)
        parfor tt3 = 1:size(finalEPIresh,2)
            multiWaitbar('Filtering & Regression', tt3/size(finalEPIresh,2));
            SignCorr = finalEPIresh(:,tt3);
            MMM = mean(SignCorr);
            
            % Detrend
            SignCorr = detrend(SignCorr,'linear');  %detrend with breakpoints over 128 seconds
            SignCorr = SignCorr + MMM;
            
            % Regression. basically includes: Movement regressors(6 parameters), mean White Matter signal and mean CSF signal
            b = regress(SignCorr,Frgval); 
            WReg  = Frgval*b;
            s_corr = SignCorr - WReg;
            s_corr = s_corr + b(end);
            
            finalEPIreshTMP(:,tt3) = s_corr;
        end
        finalEPIresh = finalEPIreshTMP;
    end
    if isequal(Regstr,'') && isequal(BinFazFilt,1)
        if ~isempty(HdLOW)
            HdLOW_Numerator = HdLOW.Numerator;
            fLPF = 1;
        else
            fLPF = 0;
        end
        if ~isempty(HdHIGH)
            HdHIGH_Numerator = HdHIGH.Numerator;
            fHPF = 1;
        else
            fHPF = 0;
        end
        parfor tt3 = 1:size(finalEPIresh,2)
            multiWaitbar('Filtering & Regression', tt3/size(finalEPIresh,2));
            SignCorr = finalEPIresh(:,tt3);
            MMM = mean(SignCorr);
            
            % Detrend
            SignCorr = detrend(SignCorr,'linear');  %detrend with breakpoints over 128 seconds
            SignCorr = SignCorr + MMM;
            
            % Low Pass Filtering
            if fLPF % Case of no low pass filter
                if NumT1<=tmp1 % if the n of dynamics if lower than 3 times the size of the filtering object
                    s_filt2 = repmat(SignCorr,ceil(tmp1/NumT1),1); % we use the center copy, avoinding phase distortion
                else
                    s_filt2 = repmat(SignCorr,3,1);
                end
                s_filtAB1LOW = filtfilt(HdLOW_Numerator,aFil,s_filt2);

                s_filtAB1LOW = s_filtAB1LOW((1+NumT1):(2*NumT1));
            else
                s_filtAB1LOW = SignCorr;
            end

            % High Pass Filtering
            if fHPF % Case of no high pass filter
                if NumT1<=tmp2
                    s_filtAB1LOW = repmat(s_filtAB1LOW,ceil(tmp2/NumT1),1);
                else
                    s_filtAB1LOW = repmat(s_filtAB1LOW,3,1);
                end
                s_filtAB1LOW_HIGH = filtfilt(HdHIGH_Numerator,aFil,s_filtAB1LOW);
                s_filtAB1LOW_HIGH = s_filtAB1LOW_HIGH((1+NumT1):(2*NumT1));
            else
                s_filtAB1LOW_HIGH = s_filtAB1LOW;
            end
            finalEPIreshTMP(:,tt3) = s_filtAB1LOW_HIGH(:) + MMM;
        end
        finalEPIresh = finalEPIreshTMP;
    end
    multiWaitbar('Filtering & Regression',1);

    clear s_filtAB1LOW s_filtAB1LOW_HIGH s_filtABF s_filt2 s_corr s_new file1 matrix4 c2Map c3Map

    finalEPIresh2 = bsxfun(@and,vetBin,ones(1,Dt)');
    finalEPIresh2 = double(finalEPIresh2);
    finalEPIresh2(finalEPIresh2~=0) = finalEPIresh;

    finalEPI = reshape(finalEPIresh2',Dx,Dy,Dz,Dt);
    finalEPI = squeeze(finalEPI);
    
    matTR = reshape(mean(finalEPI,4),prod([Dx,Dy,Dz],1));
    vcTester = any(and(0<matTR,matTR<1));
    
    if vcTester
        finalEPI = 1000.*finalEPI;
        disp('Factor 1000 applyed to the final image')
    end
    
    fprintf('Done! ============== %s\r\n',spm('time'));
end