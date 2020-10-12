%% Part of a batch process Required after CER_analysis_03
%% Analyses complex spikes in isolation
%% PRIMARILY for analysing CS and their effect on SS

% Created by NAVEEN ON 06/12/17 at CUMC

% function CER_analysis_04_LINEAR

% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------
% codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
% data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data';


cd(data_dir)
disp('!!! CER_analysis_04_ADVANCED has started running !!!');

cd(PathName1);
load(filenm);


%%
%% ONE: COMPLEX SPIKES ----------------------------------------------------


disp('*************************************************');
disp('*************************************************');
disp('!!! Analysing COMPLEX spikes !!!');



% PART 1: RASTERS --------------------

ff = figure();

suptitle('Rasters for complex spikes');

clear F;
subplot(2,2,1)
% F = Raster_n(Spikes.S,Infos(:,4),-100,1000,[0 0 0],2,Infos(:,11)); % Aligned to target
cd(codes_dir);
Raster_n(Spikes.C,Infos(:,4),-450,1250,[1 0.5 0.5],1.25); % Aligned to target
hold on;
plot(xlim,[CHANGE CHANGE],'-k','LineWidth',2);
hold off;
title('Aligned to target');

clear F;
subplot(2,2,2)
% F = Raster_n(Spikes.S,Infos(:,4),-100,1000,[0 0 0],2,Infos(:,11)); % Aligned to target
cd(codes_dir);
Raster_n(Spikes.C,Infos(:,11),-950,750,[1 0.5 0.5],1.25); % Aligned to movement
hold on;
plot(xlim,[CHANGE CHANGE],'-k','LineWidth',2);
hold off;
title('Aligned to reward');


filename = 'Rasters_Complex Spikes';
cd(Results_dir);
print(ff, '-dpdf', filename, '-r400')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: PSTHs --------------------
% All trials, @T, @R

Sigma = 100;
clear ff
ff = figure('PaperOrientation','Portrait','Units','Centimeters','paperunits','centimeters','Papertype','usletter',...
    'paperposition',[0.63452 0.64732 20.305 26.624],'Position',[5.7 1.3 15.2 17.5]); %,'Color','w');
clear F;
subplot(4,2,1)
cd(codes_dir);
P1 = HISTOGRAM_n(Spikes.C,Infos(:,4),-450,1250,Sigma,[0.2118    0.3922    0.5451]); % Aligned to target
title('Aligned to target');
% ylim([50 150])
xLim = xlim;
xlim([xLim(1)+50 xLim(2)-50])
% YLIM_1 = [min(P1)-10 max(P1)+10];
YLIM_1 = ylim;
clear ylim;
xlabel([])

subplot(4,2,2)
clear xlim xLim P;
cd(codes_dir);
P2 = HISTOGRAM_n(Spikes.C,Infos(:,11),-950,750,Sigma,[0.2353    0.7020    0.4431]); % Aligned to reward
title('Aligned to reward');
% YLIM_2 = [min(P2)-10 max(P2)+10];
YLIM_2 = ylim;
clear ylim;
xLim = xlim;
xlim([xLim(1)+50 xLim(2)-50])
xlabel([])

% AFTER LEARNT

subplot(4,2,3)
Align_time = Infos(LEARNT:size(Infos,1),4);
Start_time = -450;
End_time = 1250;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(LEARNT:size(Infos,1));
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end
% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_3 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_3 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
ylabel([])
title('After learning','fontsize',10);


subplot(4,2,4)
Align_time = Infos(LEARNT:size(Infos,1),11);
Start_time = -750;
End_time = 950;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(LEARNT:size(Infos,1));
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end
% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_4 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_4 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
ylabel([])
title('After learning','fontsize',10);


% DURING LEARNING

subplot(4,2,5)
Align_time = Infos(CHANGE:CHANGE+15,4);
Start_time = -750;
End_time = 1250;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(CHANGE:CHANGE+15);
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end
% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_5 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_5 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
ylabel([])
title('During learning','fontsize',10);

subplot(4,2,6)
Align_time = Infos(CHANGE:CHANGE+15,11);
Start_time = -750;
End_time = 950;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(CHANGE:CHANGE+15);
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end

% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_6 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_6 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
ylabel([])
title('During learning','fontsize',10);


% BEFORE CHANGE
subplot(4,2,7)
Align_time = Infos(1:CHANGE,4);
Start_time = -750;
End_time = 1250;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(1:CHANGE);
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end
% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_7 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_7 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
title('Before learning','fontsize',10);


subplot(4,2,8)
Align_time = Infos(1:CHANGE,11);
Start_time = -750;
End_time = 950;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.C(1:CHANGE);
Sig_L = SIGNAL;
Sig_R = SIGNAL;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end
% P_Left  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
% P_Right = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
% ylim_L = [min(P_Left) max(P_Left)];
% ylim_R = [min(P_Right) max(P_Right)];
% YLIM_8 = ([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))]);

P_TOTAL = HISTOGRAM_n(SIGNAL,Align_time,Start_time,End_time,Sigma,[0 0 0]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
YLIM_8 = ([min(P_TOTAL) max(P_TOTAL)]);
xlabel([])
title('Before learning','fontsize',10);


for kkk=1:8
    subplot(4,2,kkk)
    ylim([  min([YLIM_1(1) YLIM_2(1) YLIM_3(1) YLIM_4(1) YLIM_5(1) YLIM_6(1) YLIM_7(1) YLIM_8(1)]) ...
        max([YLIM_1(2) YLIM_2(2) YLIM_3(2) YLIM_4(2) YLIM_5(2) YLIM_6(2) YLIM_7(2) YLIM_8(2)])  ]);
    hold on;
    plot([0 0],ylim,'-k','linewidth',1);
end


suptitle(strcat('HIS-CS-',FileName1(6:13),'-',FileName1(15:17) ));
cd(Results_dir)
filename = strcat(FileName1(6:17),'_PSTH_Complex spike');
print(ff, '-dpdf', filename, '-r400')















