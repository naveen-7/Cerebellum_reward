%% Makes the cell plate and gives all necessary information about the CELL
%% Can be an alternative to CELL_PROPERTIES
%% BUT ONLY for display purposes
% Created by NAVEEN ON 12/05/16 at CUMC


% %
% % % function Prelim_FIRST_n
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
disp('!!! Prelim_n has started running !!!');






% Running CODE1 ------------------------
cd(PathName1);
load(filenm);

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


% FF = figure('PaperOrientation','Portrait','Units','Centimeters','paperunits','centimeters','Papertype','usletter',...
%     'paperposition',[0.63452 0.64732 20.305 26.624],'Position',[5.7 1.3 15.2 17.5]); %,'Color','w');


FF = figure();

iter=0;
ROWS=double(10);
SUB = 1;
hh = suptitle(strcat(FileNAME(1:8),'-',FileNAME(10:length(FileNAME))));
set(hh,'FontSize',10,'FontWeight','bold')







% LEARNING CURVE and RT---------------------------------------------------

n_bin = 0;  % number of bins
w_bin = round(0.1*(size(Infos,1))); % bin width
s_bin = round(0.05*(size(Infos,1)));  % shift in bin

clear per_corr1 x_trial1 RT_bef
clear per_corr2 x_trial2 RT_aft
x_trial1 = [];
x_trial2 = [];
per_corr1 = [];
per_corr2 = [];



% Part 1: before CHANGE

m_bin = round((CHANGE-1)/(s_bin))-1; % max number of bins


for i=1:s_bin:CHANGE-1
    if(i<=(m_bin-1)*s_bin)
        n_bin = n_bin+1;
        count = 0;
        RT_bef(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        RT_std_bef(n_bin) = nanstd(Infos(i:i+w_bin-1,14))/sqrt(s_bin);
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr1(n_bin) = (count/w_bin)*100;
        x_trial1(n_bin) = (i+(i+w_bin-1))/2;
    end
end


% Part 2: after CHANGE
n_bin = 0;
m_bin = round((size(Infos,1)-CHANGE)/(s_bin))-1; % max number of bins
% clear per_corr2 x_trial2

for i=CHANGE:s_bin:size(Infos,1)
    if(i<=((m_bin-2)*s_bin)+CHANGE)
        n_bin = n_bin+1;
        count = 0;
        RT_aft(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        RT_std_aft(n_bin) = nanstd(Infos(i:i+w_bin-1,14))/sqrt(s_bin);
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr2(n_bin) = (count/w_bin)*100;
        x_trial2(n_bin) = (i+(i+w_bin-1))/2;
        
        
    end
end


per_corr = NaN(length(per_corr1)+length(per_corr2),1);
per_corr(1:length(per_corr1)) = per_corr1;
per_corr(length(per_corr1)+1:length(per_corr1)+length(per_corr2)) = per_corr2;

x_trial = NaN(length(x_trial1)+length(x_trial2),1);
x_trial(1:length(x_trial1)) = x_trial1;
x_trial(length(x_trial1)+1:length(x_trial1)+length(x_trial2)) = x_trial2;

RT_tot = NaN(length(RT_bef)+length(RT_aft),1);
RT_tot(1:length(RT_bef)) = RT_bef;
RT_tot(length(RT_bef)+1:length(RT_bef)+length(RT_aft)) = RT_aft;

RT_std_tot = NaN(length(RT_std_bef)+length(RT_std_aft),1);
RT_std_tot(1:length(RT_std_bef)) = RT_std_bef;
RT_std_tot(length(RT_std_bef)+1:length(RT_std_bef)+length(RT_std_aft)) = RT_std_aft;

RT_tot_LC = RT_tot;
RT_std_tot_LC = RT_std_tot;
per_corr_LC = per_corr;
x_trial_LC = x_trial;


a = (5*SUB)+1; b = 5*(SUB+1)+1;
subplot (ROWS,5,[double(a) double(b)])

hold on;
clear ylim yLim
plot(x_trial_LC,per_corr_LC,'-','color',[1	0.7255	0.0588],'LineWidth',2);
plot(x_trial_LC,per_corr_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
xlabel('Trial number','FontSize',7);
ylabel('% Correct','FontSize',7);
set(gca,'fontsize',7)
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
plot(x_trial_LC,RT_tot_LC,'-','color',[0.8235    0.4118    0.1176],'LineWidth',2);
plot(x_trial_LC,RT_tot_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
e = errorbar(x_trial_LC,RT_tot_LC,RT_std_tot_LC);
e.Color = [0.8235    0.4118    0.1176];
xlabel('Trial number','FontSize',7);
ylabel('RT in ms','FontSize',7);
set(gca,'fontsize',7)
ylim([min(RT_tot_LC-RT_std_tot_LC) max(RT_tot_LC+RT_std_tot_LC)])
xlim([x_trial(1) x_trial(length(x_trial))]);
title('RT','FontSize',8);
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
set(gca,'FontSize',7,'LineWidth',1)
ylabel([]);
set(gca,'fontsize',7)
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
set(gca,'fontsize',7)
box off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for run=1:2;

    
    if run==1
        SUB=SUB+1;
        SIGNAL_NOW = Spikes.S;
        SHAPE_NOW = Shapes.S;
    end
    
    if run==2
        SUB=SUB+3;
        SIGNAL_NOW = Spikes.C;
        SHAPE_NOW = Shapes.C;
    end
    
    
    
    % RASTER @T----------------------------------------------------
    clear xlim;
    subplot (ROWS,5,double(5*(SUB+2)+1));
    if run==1
        title('@ target','FontSize',8,'fontweight','bold');
    end
    Raster_n(SIGNAL_NOW,Infos(:,4),-200,1200,[0.5 0.5 0.5],0.02);
    set(gca,'FontSize',7,'LineWidth',3);
    xlabel([])
    axis off;
    set(gca,'fontsize',7)
    
    
    
    % RASTER @M----------------------------------------------------
    clear xlim;
    subplot (ROWS,5,double(5*(SUB+2)+2));
    if run==1
        title('@ movement','FontSize',8,'fontweight','bold');
    end
    Raster_n(SIGNAL_NOW,Infos(:,11),-800,600,[0.5 0.5 0.5],0.02);
    set(gca,'FontSize',7,'LineWidth',3);
    set(gca,'fontweight','bold')
    xlabel([])
    axis off;
    set(gca,'fontsize',7)
    
    
    
    
    
    
    % HAND_@T -------------------------------------------------------------
    subplot (ROWS,5,double(5*(SUB+3)+1));
    Align_time = Infos(:,4);
    Start_time = -250;
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
        xlabel('time in ms','FontSize',7) 
    end
    if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
        ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
    end
    set(gca,'FontSize',7,'LineWidth',1)
    ylabel([]);
    set(gca,'fontsize',7)
    
    YLIM_T = ylim;
    clear ylim;
    
    
    % HAND_@R -------------------------------------------------------------
    subplot (ROWS,5,double(5*(SUB+3)+2));
    Align_time = Infos(:,11);
    Start_time = -850;
    End_time = 650;
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
    ylabel('Firing Rate Hz','FontSize',7)
    if run==1 
        xlabel([]) 
    end
    if run==2 
        xlabel('time in ms','FontSize',7) 
    end
    if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
        ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
    end
    set(gca,'FontSize',7,'LineWidth',1)
    ylabel([]);
    set(gca,'fontsize',7)
    
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
    subplot (ROWS,5,[double(5*(SUB+2)+3) double(5*(SUB+3)+3) ]);
    
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
        
        
        [nb,xb] = hist(TEMP_new);
        bh = bar(xb,nb);
        set(bh,'FaceColor',[0.1255    0.6980    0.6667]);
        if run==2 xlabel('ISI (in ms)','FontSize',7); end
        if ~isempty(TEMP_new)
            xlim([0 max(xb)]);
            ylim([min(nb) max(nb)]);
        end
%         if run==1 title('ISI','FontSize',8,'fontweight','bold'); end
        set(gca,'FontSize',7,'LineWidth',1)
        ylabel([]);
        set(gca,'fontsize',7)
        box off;
        
    end
    
    
    
    % SHAPE --------------------------------
    subplot (ROWS,5,[double(5*(SUB+2)+4) double(5*(SUB+3)+4) ]);
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
    
   
end





























disp(' >>> Printing the file <<<');

filename = strcat(FileNAME,'_','_DETAILS');
%  set(gcf,'PaperOrientation','Portrait','Units','Centimeters','paperunits','centimeters','Papertype','usletter',...
%     'paperposition',[0.63452 0.64732 20.305 26.624],'Position',[5.7 1.3 15.2 17.5]); %,'Color','w');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400')
%delete(gcf)




disp('!!! END OF CODE Prelims_n !!!');
