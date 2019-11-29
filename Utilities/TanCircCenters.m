function [C1,C2] = TanCircCenters(P1,P2,Ratio,plotB)
% Brunno Machado de Campos
% University of Campinas, 2017
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


% Example
% P1 = [10,10]; % Vector for point 1
% P2 = [18,18]; % Vector for point 2
% Ratio = 30; % Ratio size
% plotB = 1; % Binary for Plotting: 1 plot figure, 2 no.

Fat = [0,0,0];
if any(P1<1)
    Fat(1) = 10.*(1/P1(find(P1<1,1)));
end
if any(P2<1)
    Fat(2) = 10.*(1/P2(find(P2<1,1)));
end
if Ratio<1
    Fat(3) = 10.*(1/Ratio);
end

if any(Fat>0)
   Fat = max(Fat); 
   P1 = Fat.*P1;
   P2 = Fat.*P2;
   Ratio = Fat.*Ratio;
end

RAIO = Ratio;

ang = 0:0.002:2*pi; 
xp1 = (RAIO*cos(ang)) + P1(1);
yp1 = (RAIO*sin(ang)) + P1(2);
xp2 = (RAIO*cos(ang)) + P2(1);
yp2 = (RAIO*sin(ang)) + P2(2);
Pdis = zeros(1,size(xp1,2));
for i = 1:size(xp1,2)
    Pdis(i) = min(abs(xp1(i) - xp2) + abs(yp1(i) - yp2));
end

Pdis = round(Pdis*1000000000)./1000000000;
DPdis = diff(Pdis)./abs(diff(Pdis));
DPdis(isnan(DPdis)) = 0;
IDxs = strfind(DPdis,[-1 -1 1]);

if any(Fat>0)
    C1 = [xp1(IDxs(1)),yp1(IDxs(1))]./Fat;
    C2 = [xp1(IDxs(2)),yp1(IDxs(2))]./Fat;
else
    C1 = [xp1(IDxs(1)),yp1(IDxs(1))];
    C2 = [xp1(IDxs(2)),yp1(IDxs(2))];
end

if plotB
    if any(Fat>0)
        P1 = P1./Fat;
        P2 = P2./Fat;
        RAIO = RAIO./Fat;
        xp1 = (RAIO*cos(ang)) + P1(1);
        yp1 = (RAIO*sin(ang)) + P1(2);
        xp2 = (RAIO*cos(ang)) + P2(1);
        yp2 = (RAIO*sin(ang)) + P2(2);
    end
    
    figure
    plot(xp1,yp1 ,'LineWidth',.3,'Color',[0.9 0.9 0.9],'LineStyle','--');
    hold on
    plot(xp2,yp2 ,'LineWidth',.3,'Color',[0.9 0.9 0.9],'LineStyle','--');
    scatter(P1(1),P1(2),'MarkerEdgeColor',[0.5 0.5 0.5])
    scatter(P2(1),P2(2),'MarkerEdgeColor',[0.5 0.5 0.5])
    scatter(C1(1),C1(2),'Marker','+','MarkerEdgeColor','r','SizeData',100)
    scatter(C2(1),C2(2),'Marker','+','MarkerEdgeColor','r','SizeData',100)
    xp12 = (RAIO*cos(ang)) + C1(1);
    yp12 = (RAIO*sin(ang)) + C1(2);
    xp22 = (RAIO*cos(ang)) + C2(1);
    yp22 = (RAIO*sin(ang)) + C2(2);
    plot(xp12,yp12 ,'LineWidth',.5,'Color',[0 0 0],'LineStyle','-');
    plot(xp22,yp22 ,'LineWidth',.5,'Color',[0 0 0],'LineStyle','-');
    xlim([min([P1(1),P2(1)])-(RAIO*2) max([P1(1),P2(1)])+(RAIO*2)])
    ylim([min([P1(2),P2(2)])-(RAIO*2) max([P1(2),P2(2)])+(RAIO*2)])
    daspect([1 1 1])
    text(C1(1),C1(2),[' Center 1 (',num2str(round(C1(1)*100)/100),',',num2str(round(C1(2)*100)/100),')'])
    text(C2(1),C2(2),[' Center 2 (',num2str(round(C2(1)*100)/100),',',num2str(round(C2(2)*100)/100),')'])
end
end