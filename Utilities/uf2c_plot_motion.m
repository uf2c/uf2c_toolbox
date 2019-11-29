function [fig,HDvV,HTvV] = uf2c_plot_motion(RP_File,Visible)
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

% This function plots the SPM realigment parameters (rp_**.txt)
% 
% Input (RP_File) the whole name (path+filename) to the rp_****.txt file
% Visible = string, 'on' or 'off', for the visualization of the plot

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Ploting motion parameters\r\n');

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);
clear vetore vetoretmp
vetore = zeros(1,6);
if ischar(RP_File)
    vetore = importdata(RP_File);
else
    for yg = 1:numel(RP_File)
        vetoretmp = importdata(RP_File{yg});
        vetore = vertcat(vetore,vetoretmp(:,1:6));
        clear vetoretmp
    end
    vetore(1,:) = [];
end
vetore(:,4:6) = (180/pi).*vetore(:,4:6);

fig = figure('Visible',Visible);
set(fig,'Name','Realignment Parameters',...
    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1).*0.6 ScreSize(1).*0.4]),...
                'Color',[1 0.94 0.86]);

figTra = subplot(2,1,1);
plot(vetore(:,1))
hold on
plot(vetore(:,2))
hold on
plot(vetore(:,3))

title('Translation movement')
legend('X','Y','Z','Location','southwest')
ylabel('mm')
xlabel('volume')

figRot = subplot(2,1,2);
plot(vetore(:,4))
hold on
plot(vetore(:,5))
hold on
plot(vetore(:,6))

title('Rotational movement')
legend('Row','Pitch','Yaw','Location','southwest')
ylabel('degrees')
xlabel('volume')

MaxX = max(abs(vetore(:,1)));
MaxY = max(abs(vetore(:,2)));
MaxZ = max(abs(vetore(:,3)));

MaxR = max(abs(vetore(:,4)));
MaxP = max(abs(vetore(:,5)));
MaxYa = max(abs(vetore(:,6)));

MaxT = max([MaxX,MaxY,MaxZ]);
MaxRot = max([MaxR,MaxP,MaxYa]);

switch  MaxT
    case MaxX
        strHd = 'X';
        PosiDH = find(abs(vetore(:,1))==MaxX);
    case MaxY
        strHd = 'Y';
        PosiDH = find(abs(vetore(:,2))==MaxY);
    case MaxZ
        strHd = 'Z';
        PosiDH = find(abs(vetore(:,3))==MaxZ);
end

MaxT = num2str(MaxT);
try
    HDvV = [MaxT(1:5) ' on '   strHd];
catch
    HDvV = [MaxT ' on '   strHd];
end

switch  MaxRot
    case MaxR
        strHt = 'Row';
        PosiTH = find(abs(vetore(:,4))==MaxR);
    case MaxP
        strHt = 'Pitch';
        PosiTH = find(abs(vetore(:,5))==MaxP);
    case MaxYa
        strHt = 'Yaw';
        PosiTH = find(abs(vetore(:,6))==MaxYa);
end
MaxRot = num2str(MaxRot);

try
    HTvV = [MaxRot(1:5) ' on '   strHt];
catch
    HTvV = [MaxRot ' on '   strHt];
end



