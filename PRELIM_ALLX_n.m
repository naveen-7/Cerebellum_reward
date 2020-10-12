%% Makes the cell plate and gives all necessary information about the CELL
%% BUT ONLY for display purposes
% Created by NAVEEN ON 8/11/17 at CUMC


% %
% % % function PRELIM_ALLX_n
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
disp('!!! PRELIM_ALLX_n has started running !!!');






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


% LEARNING CURVE ---------------------------------------------------
    %     SUB=SUB-1;
    a = (5*SUB)+1; b = 5*(SUB+1)+1;
    subplot (ROWS,5,[double(a) double(b)])
    
    % Constructing the learning curve --------------
    % Method 1: continuous --------
    
    w_bin = round(0.1*(size(Infos,1))); % bin width
    s_bin = round(0.05*(size(Infos,1)));  % shift in bin
    n_bin = 0;  % number of bins
    
    m_bin = round(size(Infos,1)/(s_bin))-1; % max number of bins
    
    clear per_corr x_trial
    
    for i=1:s_bin:size(Infos,1)
        if(i<=(m_bin-1)*s_bin)
            n_bin = n_bin+1;
            count = 0;
            for j=i:i+w_bin-1
                if Infos(j,10)==1
                    count = count+1;
                end
            end
            per_corr(n_bin) = (count/w_bin)*100;
            x_trial(n_bin) = (i+(i+w_bin-1))/2;
            
        end
    end
    
    hold on;
    clear ylim yLim
    plot(x_trial,per_corr,'-','color',[1	0.7255	0.0588],'LineWidth',2);
    plot(x_trial,per_corr,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
    xlabel('Trial number','fontweight','bold','FontSize',7);
    ylabel('% Correct','fontweight','bold','FontSize',7);
    set(gca,'fontsize',7)
    
    
    yLim = ylim;
    if yLim(1)>20
        ylim([20 yLim(2)]);
    end
    yLim = ylim;
    xlim([x_trial(1) x_trial(length(x_trial))]);
    
    
    if iter==1 title('Learning curve','fontweight','bold','FontSize',8); end
    set(gca,'LineWidth',1)
    
    
    
    disp(' >>> Hold on for a bit <<<');
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run=1;


if run==1
    SUB=SUB+1;
    SIGNAL_NOW = Spikes.S;
    SHAPE_NOW = Shapes.S;
end


% RASTER @FP----------------------------------------------------
clear xlim;
s5 = subplot (ROWS,5,double(5*(SUB+2)+1));
s5Pos = get(s5,'position');
Raster_n(SIGNAL_NOW,Infos(:,3),-450,1250,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% RASTER @T----------------------------------------------------
clear xlim;
s6 = subplot (ROWS,5,double(5*(SUB+2)+2));
s6Pos = get(s6,'position');
Raster_n(SIGNAL_NOW,Infos(:,4),-450,1250,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)


% RASTER @M----------------------------------------------------
clear xlim;
s7 = subplot (ROWS,5,double(5*(SUB+2)+3));
s7Pos = get(s7,'position');
Raster_n(SIGNAL_NOW,Infos(:,11),-750,950,[0.5 0.5 0.5],0.02);
xlabel([])
axis off;
set(gca,'fontsize',5,'LineWidth',1)




% HAND_@FP -------------------------------------------------------------
s8 = subplot (ROWS,5,double(5*(SUB+3)+1));
s8Pos = get(s8,'position');
Align_time = Infos(:,3);
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
    xlabel('FP on','FontSize',7)
end
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',5,'LineWidth',0.5)
ylabel([]);

YLIM_T = ylim;
clear ylim;





% HAND_@T -------------------------------------------------------------
s9 = subplot (ROWS,5,double(5*(SUB+3)+2));
s9Pos = get(s9,'position');
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
s10 = subplot (ROWS,5,double(5*(SUB+3)+3));
s10Pos = get(s10,'position');
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


 a = 7; b = 12;
    subplot (ROWS,5,[double(a) double(b)])
    

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
 a = 8; b = 13;
    subplot (ROWS,5,[double(a) double(b)])
    

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
        subplot (ROWS,5,[double(a) double(b)])
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





set(s5,'position',s5Pos);
set(s6,'position',s6Pos);
set(s7,'position',s7Pos);
set(s8,'position',s8Pos);
set(s9,'position',s9Pos);
set(s10,'position',s10Pos);






disp(' >>> Printing the file <<<');

filename = strcat(NOME,'_','_DETAILS');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400')





disp('!!! END OF CODE PRELIM_ALLX_n !!!');
