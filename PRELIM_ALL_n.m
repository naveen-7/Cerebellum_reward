%% Makes the cell plate and gives all necessary information about the CELL
%% BUT ONLY for display purposes
% Created by NAVEEN ON 7/23/17 at CUMC


% %
% % % function PRELIM_ALL_n
% % %
% % % clc;
% % % clear all;
% % % close all;
% %
% %
% % % % Setup directories--------------------------------------------------------
% %
% % % codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW');
% % % data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');
% %
% % cd(data_dir)
disp('!!! PRELIM_ALL_n has started running !!!');






% Running CODE1 ------------------------
cd(PathName1);
load(ALLCELLS_file);

if exist('Shapes')==0
    Shapes.S = cell(size(Spikes.S));
    Shapes.C = cell(size(Spikes.C));
end

%% %%%%%%% printing the CELL PLATE %%%%%%%%% %%

% CELL INFORMATION---------------------------------------------

clear CC matches;
str = FileName1(6:length(FileName1)-4);
[CC,matches] = strsplit(str,{'_'});
FileNAME = str;


FF = figure();

iter=0;
ROWS=double(10);
SUB = 1;
hh = suptitle(strcat(NOME(1:8),'-',NOME(10:length(NOME))));
set(hh,'FontSize',10,'FontWeight','bold')


% LEARNING CURVE and RT---------------------------------------------------

W_bin = 0.1; % bin width
S_bin = 0.05;  % shift in bin
[x_trial, RT_trial, per_corr] = LEARNING_CURVE_n(Infos(:,14), CHANGE, Infos, W_bin, S_bin, 'p');

a = (5*SUB)+1; b = 5*(SUB+1)+1;
subplot (ROWS,5,[double(a) double(b)])

hold on;
clear ylim yLim
plot(x_trial,per_corr,'-','color',[0.5 0.5 0.5],'LineWidth',1.5);
plot(x_trial,per_corr,'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
xlabel('Trial number','FontSize',7);
ylabel('% Correct','FontSize',7);
set(gca,'fontsize',5)
yLim = ylim;
if yLim(1)>20
    ylim([20 yLim(2)]);
end
yLim = ylim;
xlim([x_trial(1) x_trial(length(x_trial))]);
title('Learning Curve','FontSize',8);
set(gca,'LineWidth',1)
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);

a = (5*SUB)+2; b = 5*(SUB+1)+2;
subplot (ROWS,5,[double(a) double(b)])
hold on;
clear ylim yLim
plot(x_trial,RT_trial(:,1),'-','color',[1 0.4 0.4],'LineWidth',1.5);
errorbar(x_trial,RT_trial(:,1),RT_trial(:,3),'.k');
plot(x_trial,RT_trial(:,1),'o','MarkerSize',3,'MarkerFaceColor',[1 0.4 0.4],'color','k');
xlabel('Trial number','FontSize',8);
% ylabel('RT','FontSize',8);
title('RT','FontSize',8);
set(gca,'fontsize',5)

ylim([min(RT_trial(:,1)-RT_trial(:,3)) max(RT_trial(:,1)+RT_trial(:,3))])
xlim([x_trial(1) x_trial(length(x_trial))]);
set(gca,'LineWidth',1)
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);



% CS-SS ISI

% Raster
a = (5*SUB)+4; b = 5*(SUB+1)+4;
subplot (ROWS,5,[double(a) double(b)])

Start=-50;
End = 100;
CER_CS_RASTER_ALIGNED(Spikes.S,Spikes.C,Start,End)

xlabel('time from CS (in ms)','FontSize',7);
title('CS triggered SS','FontSize',8,'fontweight','bold');
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);
box off;


% distribution of ISI between CS and SS
a = (5*SUB)+3; b = 5*(SUB+1)+3;
subplot (ROWS,5,[double(a) double(b)])
Start=-10;
End = 100;
[t ISI_CS_SS] = Align_CS_n(Spikes.S,Spikes.C,Start,End);


[nb,xb] = hist(ISI_CS_SS,10);
bh = bar(xb,nb);
set(bh,'FaceColor',[0.1255    0.6980    0.6667]);
% xlabel('ISI (in ms)','FontSize',7);
if ~isnan(ISI_CS_SS)
    xlim([0 max(xb)]);
    ylim([min(nb) max(nb)]);
end
title('CS-SS-ISI','FontSize',8,'fontweight','bold');
set(gca,'FontSize',7,'LineWidth',1)
ylabel([]);
set(gca,'fontsize',5)
box off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run=1;


if run==1
    SUB=SUB+1;
    SIGNAL_NOW = Spikes.S;
    SHAPE_NOW = Shapes.S;
end

%     if run==2
%         SUB=SUB+4;
%         SIGNAL_NOW = Spikes.C;
%         SHAPE_NOW = Shapes.C;
%     end



% RASTER @T----------------------------------------------------
clear xlim;
s5 = subplot (ROWS,5,double(5*(SUB+2)+1));
s5Pos = get(s5,'position');
Raster_n(SIGNAL_NOW,Infos(:,4),-450,1250,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% RASTER @M----------------------------------------------------
clear xlim;
s6 = subplot (ROWS,5,double(5*(SUB+2)+2));
s6Pos = get(s6,'position');
Raster_n(SIGNAL_NOW,Infos(:,11),-750,950,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% HAND_@T -------------------------------------------------------------
s7 = subplot (ROWS,5,double(5*(SUB+3)+1));
s7Pos = get(s7,'position');
Align_time = Infos(:,4);
Start_time = -450;
End_time = 1250;
Sigma = 20;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
LW = 1.25;
clear SIGNAL xlim xLIM;
SIGNAL = SIGNAL_NOW;
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
P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW,0);
P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW,0);
ylim_L = [min(P_Left(50:length(P_Left-50))) max(P_Left(50:length(P_Left-500)))];
ylim_R = [min(P_Right(50:length(P_Right-50))) max(P_Right(50:length(P_Right-500)))];
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
ylabel('Firing Rate Hz','FontSize',7)
if run==1
    xlabel([])
end
if run==2
    xlabel('targ','FontSize',7)
end
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);

YLIM_T = ylim;
clear ylim;


% HAND_@R -------------------------------------------------------------
s8 = subplot (ROWS,5,double(5*(SUB+3)+2));
s8Pos = get(s8,'position');
Align_time = Infos(:,11);
Start_time = -750;
End_time = 950;
Sigma = 20;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = SIGNAL_NOW;
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
P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW,0);
P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW,0);
ylim_L = [min(P_Left(50:length(P_Left-50))) max(P_Left(50:length(P_Left-50)))];
ylim_R = [min(P_Right(50:length(P_Right-50))) max(P_Right(50:length(P_Right-50)))];
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
%     ylabel('Firing Rate Hz','FontSize',7)
if run==1
    xlabel([])
end
if run==2
    xlabel('movt','FontSize',7)
end
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);

YLIM_M = ylim;



subplot (ROWS,5,double(5*(SUB+3)+1));
ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
subplot (ROWS,5,double(5*(SUB+3)+2));
ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);



% % ISI ------------------------------------------------------
s9 = subplot (ROWS,5,[double(5*(SUB+2)+3) double(5*(SUB+3)+3) ]);
s9Pos = get(s9,'position');

clear temp TEMP TEMP_new
cellflag=0;
for i=1:size(SIGNAL_NOW,1)
    Signal_S = SIGNAL_NOW{i,1};
    Time_F = Infos(i,3);   % FP_on
    Time_S = Infos(i,4);   % Stimulus
    
    Length = Time_S-Time_F;
    
    Sig_S = Signal_S(find(Time_F<Signal_S & Signal_S<Time_S));
    Sig_S = Sig_S-Time_F;
    cell_len = length(Spikes.S{i,1});
    if cell_len>=2
        temp{i,1} = diff(Sig_S);
        cellflag=cellflag+1;
    end
end



if cellflag>0
    
    TEMP = sort(cell2mat(temp));
    temp_avg = nanmean(TEMP);
    temp_std = nanstd(TEMP);
    TEMP_new = TEMP(1:find(TEMP>=temp_avg+temp_std,1));
    
    if size(TEMP_new,1)>=15
        h = histfit(TEMP_new,[],'gamma');
        h(1).FaceColor = [0.1255    0.6980    0.6667];
        h(2).Color = [244, 75, 66]/255;
        set(gca,'FontSize',5,'LineWidth',0.5)
        ylabel([]);
        box off;
        [phat,pci] = gamfit(TEMP_new);
        k = phat(1);
        theta = phat(2);
        XLIM = xlim; YLIM = ylim;
        text(XLIM(2)*0.6,YLIM(2)*0.7,strcat('k = ',num2str(k)),'FontSize',5)
    else
        [nb,xb] = hist(TEMP_new);
        bh = bar(xb,nb);
        set(bh,'FaceColor',[0.1255    0.6980    0.6667]);
        if run==2 xlabel('ISI (in ms)','FontSize',7); end
        if ~isempty(TEMP_new)
            xlim([0 max(xb)]);
            ylim([min(nb) max(nb)]);
        end
        set(gca,'FontSize',5,'LineWidth',0.5)
        ylabel([]);
        box off;
    end
    
    
end



% SHAPE --------------------------------
s10 = subplot (ROWS,5,[double(5*(SUB+2)+4) double(5*(SUB+3)+4) ]);
s10Pos = get(s10,'position');
clear SHAPE_MAT AVG_SHAPE STD_SHAPE

SHAPE_MAT = cell2mat(SHAPE_NOW);
if ~isempty(SHAPE_MAT)
    AVG_SHAPE = nanmean(SHAPE_MAT);
    STD_SHAPE = nanstd(SHAPE_MAT);
    hold on;
    plot(AVG_SHAPE,'-','color',[0.2118    0.3922    0.5451],'linewidth',1.5);
    clear x y1 y2 X Y
    
    x = 1:size(SHAPE_MAT,2);
    X=[x,fliplr(x)];
    
    % error
    y1 =  AVG_SHAPE + 1*( STD_SHAPE );
    y2 =  AVG_SHAPE - 1*( STD_SHAPE );
    Y=[y1,fliplr(y2)];
    fi = fill(X,Y,[0.2118    0.3922    0.5451],'linestyle','none');
    set(fi,'FaceAlpha',0.2);
    
    xlim([1 size(SHAPE_MAT,2)]);
    ylim([min(min(SHAPE_MAT)) max(max(SHAPE_MAT))]);
    axis off;
    if run==1
        title('simple spikes','FontSize',8,'fontweight','bold')
    end
    if run==2
        title('complex spikes','FontSize',8,'fontweight','bold')
    end
    hold on;
end


if ~isempty(SHAPE_MAT)
    
    Average = nanmean(AVG_SHAPE);
    STDev   = nanstd(AVG_SHAPE);
    
    UPPER = Average + 0.3*STDev;
    LOWER = Average - 0.3*STDev;
    
    DataInv = 1.01*max(AVG_SHAPE) - AVG_SHAPE;
    Top_H = nanmean(AVG_SHAPE)+0.75*nanstd(AVG_SHAPE);
    Bot_H = nanmean(DataInv)+0.75*nanstd(DataInv);
    
    [Maxima,MaxIdx] = findpeaks(AVG_SHAPE,'MinPeakHeight',Top_H); %findpeaks(AVG_SHAPE);     % Maxima value and position
    [Minima1,MinIdx] = findpeaks(DataInv,'MinPeakHeight',Bot_H);      % Minima position
    Minima = AVG_SHAPE(MinIdx);                % Minima value
    
    Point1 = find(LOWER>=AVG_SHAPE | AVG_SHAPE>=UPPER,1);
    Temp_Point = max(max(MaxIdx),max(MinIdx));
    PointX = find(LOWER<=AVG_SHAPE(Temp_Point:length(AVG_SHAPE)) & AVG_SHAPE(Temp_Point:length(AVG_SHAPE))<=UPPER,1) + Temp_Point-1;
    
    if isempty(Point1)
        Point1 = 1;
    end
    
    if isempty(PointX)
        PointX = length(AVG_SHAPE);
    end
    
    SPK_duration = PointX-Point1;
    
    iPk = [Point1 PointX];   % POINTS TO DISPLAY
    subplot (ROWS,5,[double(5*(SUB+2)+4) double(5*(SUB+3)+4) ]);
    plot(iPk,AVG_SHAPE(iPk)+0.2*(max(AVG_SHAPE)-min(AVG_SHAPE)),'VK','MarkerSize',4,'MarkerFaceColor','R');
    
    if exist ('SPK_duration')
        SPK_duration = SPK_duration*0.02; %in ms
    end
    
    xLIM = xlim;
    yLIM = ylim;
    text(((xLIM(1)+xLIM(2))/2),yLIM(1)-0.1,strcat(num2str(SPK_duration),' ms'),'Color','red','FontSize',8);
    
    
end

disp(' >>> Almost DONE <<<');


% end





%% DO THE DELTA ----------------------------------------------



Type8 = NaN(length(Infos),1);

for i=2:length(Infos)
    if (Infos(i-1,10)==1)
        Type8(i,1)=1;
    end
end

IND = CHANGE-15:CHANGE;

clear PSTH_bef PSTH_aft PSTH_bef_std PSTH_aft_std
clear points_bef_T points_aft_T points_bef_R points_aft_R

SIGMA = 30;


for XX =1:2   % Running twice; each time align differently
    
    if XX ==1 Align_code = 4; end
    if XX ==2 Align_code = 11; end
    
    if Align_code == 11
        Start = -750; End = 950;
        time_M = Start:End;
    end
    
    if Align_code == 4
        Start = -450; End = 1250;
        time_T = Start:End;
    end
    
    x = 6;
    
    time = Start:End;
    
    % BEFORE CHANGE ----------
    
    clear temp_trials
    temp_trials = find(Type8(1:CHANGE,1)>=1);
    no_trials = length(temp_trials)
    no_trials_BC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,SIGMA,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,SIGMA,[0.6353    0.8039    0.3529]);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    s11 = subplot(ROWS,5,double(5*(x+2))+XX);
    s11Pos = get(s11,'position');
    hold on;
    plot(time,nanmean(PSTH_type1),'color',[0.4000    0.8039         0],'LineWidth',1);
        plot(time,nanmean(PSTH_type2),'color',[0.2745    0.5098    0.7059],'LineWidth',1);
%     PSTH_n(Spikes.S(IND,:),Infos(IND,Align_code),Start,End,SIGMA,[0.4 0.4 0.4],1,0); % Plotting the main task
    ylabel([]); xlabel([]);
    set(gca,'FontSize',5,'LineWidth',0.5)
    xlim([Start+50 End-50]);
    
    ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))+10]);
    hold on;
    plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
    
    
    % DURING CHANGE ----------
    clear temp_trials
    temp_trials = find(Type8(CHANGE:CHANGE+15,1)>=1)+CHANGE-1;
    no_trials = length(temp_trials)
    no_trials_DC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,SIGMA,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,SIGMA,[0.6353    0.8039    0.3529]);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    s12 = subplot(ROWS,5,double(5*(x+1))+XX);
    s12Pos = get(s12,'position');
    hold on;
    plot(time,nanmean(PSTH_type1),'color',[0.4000    0.8039         0],'LineWidth',1);
        plot(time,nanmean(PSTH_type2),'color',[0.2745    0.5098    0.7059],'LineWidth',1);
%     PSTH_n(Spikes.S(IND,:),Infos(IND,Align_code),Start,End,SIGMA,[0.4 0.4 0.4],1,0); % Plotting the main task
    ylabel([]); xlabel([]);
    set(gca,'FontSize',5,'LineWidth',0.5)
    xlim([Start+50 End-50]);
    
    ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))+10]);
    hold on;
    plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
    
    
    % AFTER CHANGE ----------
    
    clear temp_trials
    temp_trials = find(Type8(LEARNT+5:size(Infos,1),1)>=1)+LEARNT;
    no_trials = length(temp_trials)
    no_trials_AC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,(End-Start));
    PSTH_type2 =  NaN(no_trials,(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,SIGMA,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,SIGMA,[0.6353    0.8039    0.3529]);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    s13 = subplot(ROWS,5,double(5*(x))+XX);
    s13Pos = get(s13,'position');
    hold on;
    plot(time,nanmean(PSTH_type1),'color',[0.4000    0.8039         0],'LineWidth',1);
        plot(time,nanmean(PSTH_type2),'color',[0.2745    0.5098    0.7059],'LineWidth',1);
%     PSTH_n(Spikes.S(IND,:),Infos(IND,Align_code),Start,End,SIGMA,[0.4 0.4 0.4],1,0); % Plotting the main task
    ylabel([]); xlabel([]);
    set(gca,'FontSize',5,'LineWidth',0.5)
    xlim([Start+50 End-50]);
    
    
    ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))+10]);
    hold on;
    plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
    
    
    
    set(s11,'position',s11Pos);
    set(s12,'position',s12Pos);
    set(s13,'position',s13Pos);
    
end



set(s11,'position',s11Pos);
set(s12,'position',s12Pos);
set(s13,'position',s13Pos);


%% CS ------------------------------------------------

SIGNAL_NOW = Spikes.C;
SHAPE_NOW = Shapes.C;


% RASTER @T----------------------------------------------------
clear xlim;
s14 = subplot (ROWS,5,33);
s14Pos = get(s14,'position');
Raster_n(SIGNAL_NOW,Infos(:,4),-450,1250,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% RASTER @M----------------------------------------------------
clear xlim;
s15 = subplot (ROWS,5,34);
s15Pos = get(s15,'position');
Raster_n(SIGNAL_NOW,Infos(:,11),-750,950,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% HAND_@T -------------------------------------------------------------
s16 = subplot (ROWS,5,38);
s16Pos = get(s16,'position');
Align_time = Infos(:,4);
Start_time = -450;
End_time = 1250;
Sigma = 20;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
LW = 1;
clear SIGNAL xlim xLIM;
SIGNAL = SIGNAL_NOW;
Sig_L = SIGNAL;
Sig_R = SIGNAL;

Bin = 200;

for i=1:length(SIGNAL)
    if (Infos(i,9)==1)      % Left
        Sig_L{i,1}=[];
    end
    if (Infos(i,9)==0)      % Right
        Sig_R{i,1}=[];
    end
end

[x, P_Left]  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Bin,Colour1,LW,1);
[x, P_Right] = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Bin,Colour2,LW,1);
% P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW,0);
% P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW,0);
ylim_L = [min(P_Left) max(P_Left)];
ylim_R = [min(P_Right) max(P_Right)];
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
ylabel('Firing Rate Hz','FontSize',7)
if run==1
    xlabel([])
end
if run==2
    xlabel('targ','FontSize',7)
end
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);

YLIM_T = ylim;
clear ylim;


% HAND_@R -------------------------------------------------------------
s17 = subplot (ROWS,5,39);
s17Pos = get(s17,'position');
Align_time = Infos(:,11);
Start_time = -750;
End_time = 950;
Sigma = 20;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = SIGNAL_NOW;
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

[x, P_Left]  = HISTOGRAM_n(Sig_L,Align_time,Start_time,End_time,Bin,Colour1,LW,1);
[x, P_Right] = HISTOGRAM_n(Sig_R,Align_time,Start_time,End_time,Bin,Colour2,LW,1);

% P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW,0);
% P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW,0);
ylim_L = [min(P_Left) max(P_Left)];
ylim_R = [min(P_Right) max(P_Right)];
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
%     ylabel('Firing Rate Hz','FontSize',7)
if run==1
    xlabel([])
end
if run==2
    xlabel('movt','FontSize',7)
end
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);

YLIM_M = ylim;



subplot (ROWS,5,38);
ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
subplot (ROWS,5,39);
ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);





% % ISI ------------------------------------------------------
s18 = subplot (ROWS,5,[43]);
s18Pos = get(s18,'position');
clear temp TEMP TEMP_new
cellflag=0;
for i=1:size(SIGNAL_NOW,1)
    Signal_S = SIGNAL_NOW{i,1};
    Time_F = Infos(i,3);   % FP_on
    Time_S = Infos(i,4);   % Stimulus
    
    Length = Time_S-Time_F;
    
    Sig_S = Signal_S(find(Time_F<Signal_S & Signal_S<Time_S));
    Sig_S = Sig_S-Time_F;
    cell_len = length(Spikes.S{i,1});
    if cell_len>=2
        temp{i,1} = diff(Sig_S);
        cellflag=cellflag+1;
    end
end



if cellflag>0
    
    TEMP = sort(cell2mat(temp));
    temp_avg = nanmean(TEMP);
    temp_std = nanstd(TEMP);
    TEMP_new = TEMP(1:find(TEMP>=temp_avg+temp_std,1));
    
    if size(TEMP_new,1)>=15
        h = histfit(TEMP_new,[],'gamma');
        h(1).FaceColor = [0.1255    0.6980    0.6667];
        h(2).Color = [244, 75, 66]/255;
        set(gca,'FontSize',5,'LineWidth',0.5)
        ylabel([]);
        box off;
        [phat,pci] = gamfit(TEMP_new);
        k = phat(1);
        theta = phat(2);
        XLIM = xlim; YLIM = ylim;
        text(XLIM(2)*0.6,YLIM(2)*0.7,strcat('k = ',num2str(k)),'FontSize',5)
    else
        [nb,xb] = hist(TEMP_new);
        bh = bar(xb,nb);
        set(bh,'FaceColor',[0.1255    0.6980    0.6667]);
        if run==2 xlabel('ISI (in ms)','FontSize',7); end
        if ~isempty(TEMP_new)
            xlim([0 max(xb)]);
            ylim([min(nb) max(nb)]);
        end
        set(gca,'FontSize',5,'LineWidth',0.5)
        ylabel([]);
        box off;
    end
    
end







% SHAPE --------------------------------
s19 = subplot (ROWS,5,44);
s19Pos = get(s19,'position');
clear SHAPE_MAT AVG_SHAPE STD_SHAPE



SHAPE_MAT = cell2mat(SHAPE_NOW);
if ~isempty(SHAPE_MAT)
    AVG_SHAPE = nanmean(SHAPE_MAT);
    STD_SHAPE = nanstd(SHAPE_MAT);
    hold on;
    plot(AVG_SHAPE,'-','color',[0.2118    0.3922    0.5451],'linewidth',1.5);
    clear x y1 y2 X Y
    
    x = 1:size(SHAPE_MAT,2);
    X=[x,fliplr(x)];
    
    % error
    y1 =  AVG_SHAPE + 1*( STD_SHAPE );
    y2 =  AVG_SHAPE - 1*( STD_SHAPE );
    Y=[y1,fliplr(y2)];
    fi = fill(X,Y,[0.2118    0.3922    0.5451],'linestyle','none');
    set(fi,'FaceAlpha',0.2);
    
    xlim([1 size(SHAPE_MAT,2)]);
    ylim([min(min(SHAPE_MAT)) max(max(SHAPE_MAT))]);
    axis off;
    %         title('complex spikes','FontSize',8,'fontweight','bold')
    hold on;
end


if ~isempty(SHAPE_MAT)
    
    Average = nanmean(AVG_SHAPE);
    STDev   = nanstd(AVG_SHAPE);
    
    UPPER = Average + 0.3*STDev;
    LOWER = Average - 0.3*STDev;
    
    DataInv = 1.01*max(AVG_SHAPE) - AVG_SHAPE;
    Top_H = nanmean(AVG_SHAPE)+0.75*nanstd(AVG_SHAPE);
    Bot_H = nanmean(DataInv)+0.75*nanstd(DataInv);
    
    [Maxima,MaxIdx] = findpeaks(AVG_SHAPE,'MinPeakHeight',Top_H); %findpeaks(AVG_SHAPE);     % Maxima value and position
    [Minima1,MinIdx] = findpeaks(DataInv,'MinPeakHeight',Bot_H);      % Minima position
    Minima = AVG_SHAPE(MinIdx);                % Minima value
    
    Point1 = find(LOWER>=AVG_SHAPE | AVG_SHAPE>=UPPER,1);
    Temp_Point = max(max(MaxIdx),max(MinIdx));
    PointX = find(LOWER<=AVG_SHAPE(Temp_Point:length(AVG_SHAPE)) & AVG_SHAPE(Temp_Point:length(AVG_SHAPE))<=UPPER,1) + Temp_Point-1;
    
    if isempty(Point1)
        Point1 = 1;
    end
    
    if isempty(PointX)
        PointX = length(AVG_SHAPE);
    end
    
    SPK_duration = PointX-Point1;
    
    iPk = [Point1 PointX];   % POINTS TO DISPLAY
    subplot (ROWS,5,[44 ]);
    plot(iPk,AVG_SHAPE(iPk)+0.2*(max(AVG_SHAPE)-min(AVG_SHAPE)),'VK','MarkerSize',4,'MarkerFaceColor','R');
    
    if exist ('SPK_duration')
        SPK_duration = SPK_duration*0.02; %in ms
    end
    
    xLIM = xlim;
    yLIM = ylim;
    text(((xLIM(1)+xLIM(2))/2),yLIM(1)-0.1,strcat(num2str(SPK_duration),' ms'),'Color','red','FontSize',8);
    
    
end














set(s5,'position',s5Pos);
set(s6,'position',s6Pos);
set(s7,'position',s7Pos);
set(s8,'position',s8Pos);
set(s9,'position',s9Pos);
set(s11,'position',s11Pos);
set(s12,'position',s12Pos);
set(s13,'position',s13Pos);
set(s14,'position',s14Pos);
set(s16,'position',s16Pos);
set(s17,'position',s17Pos);
set(s18,'position',s18Pos);
set(s19,'position',s19Pos);


disp(' >>> Printing the file <<<');

filename = strcat(NOME,'_','_DETAILS');
%  set(gcf,'PaperOrientation','Portrait','Units','Centimeters','paperunits','centimeters','Papertype','usletter',...
%     'paperposition',[0.63452 0.64732 20.305 26.624],'Position',[5.7 1.3 15.2 17.5]); %,'Color','w');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400')





disp('!!! END OF CODE PRELIMS_n !!!');
