%% Code for analysing Manipulandam change


% Created by NAVEEN ON 08/24/17 at CUMC


% function CER_analysis_B_LINEAR

% % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FOR REGULAR, RELEASE HERE -----------------------------------------

clc;
clear 
close all;


codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');

cd(data_dir)
disp('!!! CER_analysis_X_LINEAR has started running !!!');



% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE "M" DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile);
cd(PathName1);






% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % FOR BATCH, RELEASE HERE -----------------------------------------
% % % 
% % % close all;
% % % File = DATAfile;
% % % FileName = Tempy;
% % % 
% % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%
%% ONE: SIMPLE SPIKES ----------------------------------------------------


disp('*************************************************');
disp('*************************************************');
disp('!!! Analysing simple spikes !!!');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: PSTHs --------------------


B_COLOUR = [0.5 0.5 0.5];
D_COLOUR = [158 125 175]/255;



F = figure();


% bars -----------------------------------------------------
s1 = subplot(5,5,1);
s1Pos = get(s1,'position');

clear SHAPE_MAT AVG_SHAPE STD_SHAPE
SHAPE_NOW = Shapes.B;  



SHAPE_MAT = cell2mat(SHAPE_NOW);
if ~isempty(SHAPE_MAT)
    
    for ii=1:size(SHAPE_MAT)
         hold on;
       plot(SHAPE_MAT(ii,:),'color',[0.7 0.7 0.7],'linewidth',0.5)  
    end
     hold on;
    AVG_SHAPE = nanmean(SHAPE_MAT);
    STD_SHAPE = nanstd(SHAPE_MAT);
    hold on;
    plot(AVG_SHAPE,'-','color',B_COLOUR,'linewidth',1.2);
    clear x y1 y2 X Y
 
    xlim([1 size(SHAPE_MAT,2)]);
    ylim([min(min(SHAPE_MAT)) max(max(SHAPE_MAT))]);
    axis off;
   
    hold on;
end






% dowels -----------------------------------------------------
s2 = subplot(5,5,6);
s2Pos = get(s2,'position');

clear SHAPE_MAT AVG_SHAPE STD_SHAPE
SHAPE_NOW = Shapes.D;  



SHAPE_MAT = cell2mat(SHAPE_NOW);
if ~isempty(SHAPE_MAT)
    
    for ii=1:size(SHAPE_MAT)
        hold on;
       plot(SHAPE_MAT(ii,:),'color',[0.8 0.8 0.8],'linewidth',0.5) 
    end
    hold on;
    AVG_SHAPE = nanmean(SHAPE_MAT);
    STD_SHAPE = nanstd(SHAPE_MAT);
    hold on;
    plot(AVG_SHAPE,'-','color',D_COLOUR,'linewidth',1.2);
    clear x y1 y2 X Y
 
    xlim([1 size(SHAPE_MAT,2)]);
    ylim([min(min(SHAPE_MAT)) max(max(SHAPE_MAT))]);
    axis off;
   
    hold on;
end




Start_T = -750;
End_T  =1250;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =1250;
StartLIM_M = -400; EndLIM_M = 800;




subplot(5,5,[2 3])
hold on;
RasterBAR_n(Spikes.B(4:4+5),INFOS.B(4:4+5,4),Start_T,End_T,B_COLOUR,0.02);
axis off;
ylim([-1 7])
xlim([StartLIM_T EndLIM_T])

subplot(5,5,[7 8])
hold on;
RasterBAR_n(Spikes.D(4:4+5),INFOS.D(4:4+5,4),Start_T,End_T,D_COLOUR,0.02);
axis off;
ylim([-1 7])
axis off;
xlim([StartLIM_T EndLIM_T])



subplot(5,5,[4 5])
hold on;
RasterBAR_n(Spikes.B(4:4+5),INFOS.B(4:4+5,11),Start_M,End_M,B_COLOUR,0.02);
axis off;
ylim([-1 7])
xlim([StartLIM_M EndLIM_M])

subplot(5,5,[9 10])
hold on;
RasterBAR_n(Spikes.D(4:4+5),INFOS.D(4:4+5,11),Start_M,End_M,D_COLOUR,0.02);
axis off;
ylim([-1 7])
axis off;
xlim([StartLIM_M EndLIM_M])








temp_time = StartLIM_T:EndLIM_T;

subplot(5,5,[12 13 17 18])
hold on;
P1 = PSTHe_n(Spikes.B,INFOS.B(:,4),Start_T,End_T,30,B_COLOUR,1,1); 
P2 = PSTHe_n(Spikes.D,INFOS.D(:,4),Start_T,End_T,30,B_COLOUR,1,1);
temp_ind = find(StartLIM_T<=temp_time & temp_time<=EndLIM_T);
errorline_n(StartLIM_T:EndLIM_T,P1(1,temp_ind),P1(3,temp_ind),1,B_COLOUR,0.3,0,1)
errorline_n(StartLIM_T:EndLIM_T,P2(1,temp_ind),P2(3,temp_ind),1,D_COLOUR,0.3,0,1)
xlim([StartLIM_T EndLIM_T])
hold on;

MIN = nanmin([nanmin(P1(1,find(temp_time>=-1100 & temp_time<=500))-P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmin(P2(1,find(temp_time>=-1100 & temp_time<=500))-P2(3,find(temp_time>=-1100 & temp_time<=500)))])-2;

MAX = nanmax([nanmax(P1(1,find(temp_time>=-1100 & temp_time<=500))+P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmax(P2(1,find(temp_time>=-1100 & temp_time<=500))+P2(3,find(temp_time>=-1100 & temp_time<=500)))])+2;
MAX = 30;
ylim([MIN MAX])

plot([StartLIM_T StartLIM_T],[MIN MIN+20],'-k','linewidth',1)
plot([StartLIM_T StartLIM_T+200],[MIN MIN],'-k','linewidth',1)
axis off;

plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)







%%%%%%%%%%%%%%%%%%%%%%%%% movt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_time = StartLIM_M:EndLIM_M;

subplot(5,5,[14 15 19 20])
hold on;
P1 = PSTHe_n(Spikes.B,INFOS.B(:,11),Start_M,End_M,30,B_COLOUR,1,0); 
P2 = PSTHe_n(Spikes.D,INFOS.D(:,11),Start_M,End_M,30,D_COLOUR,1,0);
% temp_ind = find(StartLIM_M<=temp_time & temp_time<=EndLIM_M);
% errorline_n(StartLIM_M:EndLIM_M,P1(1,temp_ind),P1(3,temp_ind),1,B_COLOUR,0.3,0,1)
% errorline_n(StartLIM_M:EndLIM_M,P2(1,temp_ind),P2(3,temp_ind),1,D_COLOUR,0.3,0,1)
xlim([StartLIM_M EndLIM_M])
hold on;

MIN = nanmin([nanmin(P1(1,find(temp_time>=-1100 & temp_time<=500))-P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmin(P2(1,find(temp_time>=-1100 & temp_time<=500))-P2(3,find(temp_time>=-1100 & temp_time<=500)))])-2;

MAX = nanmax([nanmax(P1(1,find(temp_time>=-1100 & temp_time<=500))+P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmax(P2(1,find(temp_time>=-1100 & temp_time<=500))+P2(3,find(temp_time>=-1100 & temp_time<=500)))])+2;
MAX = 30;
ylim([MIN MAX])



plot([StartLIM_M StartLIM_M],[MIN MIN+20],'-k','linewidth',1)
plot([StartLIM_M StartLIM_M+200],[MIN MIN],'-k','linewidth',1)
axis off;

plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)





















suptitle(strcat('BD-',NOME(1:8),'-',NOME(10:12)));


cd(Results_dir)
filename = 'bars_dowels';
print(F, '-dpdf', filename, '-r400')




%% COLOUR PSTH

% Create a structure
% bars then dowels

clear spikes;
spikes.S(1:length(Spikes.B),1) = Spikes.B;
spikes.S(length(Spikes.B)+1:length(Spikes.B)+length(Spikes.D),1) = Spikes.D;

CHANGE=11;

clear infos
infos(1:size(INFOS.B,1),:) = INFOS.B;
infos(size(INFOS.B,1)+1:size(INFOS.B,1)+size(INFOS.D,1),:) = INFOS.D;


%% PSTH
% % % % % % % 
% % % % % % % for I=1:2
% % % % % % %       
% % % % % % %     if I==1 Align_code = 4; end
% % % % % % %     if I==2 Align_code = 11; end
% % % % % % %     
% % % % % % %     if Align_code == 11
% % % % % % %         Start = -950; End = 750;
% % % % % % %     end
% % % % % % %     
% % % % % % %     if Align_code == 4
% % % % % % %         Start = -450; End = 1250;
% % % % % % %     end
% % % % % % %     
% % % % % % %     sig_len = End-Start;
% % % % % % %     
% % % % % % %     
% % % % % % %     
% % % % % % %     
% % % % % % %     
% % % % % % %     %Method 2: dis-continuous --------
% % % % % % %     n_bin = 0;  % number of bins
% % % % % % %     w_bin = round(0.1*(size(infos,1))); % bin width
% % % % % % %     s_bin = round(0.02*(size(infos,1)));  % shift in bin
% % % % % % %     
% % % % % % %     
% % % % % % %     % Part 1: before CHANGE
% % % % % % %     
% % % % % % %     m_bin = round((CHANGE-1)/(s_bin))-1; % max number of bins
% % % % % % %     clear PSTH_ALL1 count
% % % % % % %     
% % % % % % %     for i=1:s_bin:CHANGE-1
% % % % % % %         if(i <= CHANGE-1+1-w_bin )
% % % % % % %             n_bin = n_bin+1;
% % % % % % %             clear Signal;
% % % % % % %             Signal = spikes.S(i:i+w_bin-1);
% % % % % % %             PSTH_ALL1(n_bin,:) = PSTH_n(Signal,infos(i:i+w_bin-1,Align_code),Start,End,10,[0.2118    0.3922    0.5451]);
% % % % % % %         end
% % % % % % %     end
% % % % % % %     
% % % % % % %     % Part 2: after CHANGE
% % % % % % %     n_bin = 0;
% % % % % % %     m_bin = round((size(infos,1)-CHANGE)/(s_bin)); % max number of bins
% % % % % % %     clear PSTH_ALL2 count
% % % % % % %     
% % % % % % %     for i=CHANGE:s_bin:size(infos,1)
% % % % % % %         if(i <= size(infos,1)+1-w_bin )                  %((m_bin-2)*s_bin+CHANGE))
% % % % % % %             n_bin = n_bin+1;
% % % % % % %             clear Signal;
% % % % % % %             Signal = spikes.S(i:i+w_bin-1);
% % % % % % %             PSTH_ALL2(n_bin,:) = PSTH_n(Signal,infos(i:i+w_bin-1,Align_code),Start,End,10,[0.2118    0.3922    0.5451]);
% % % % % % %         end
% % % % % % %     end
% % % % % % %     
% % % % % % %    
% % % % % % %     
% % % % % % %     clear PSTH_ALL
% % % % % % %     PSTH_ALL = NaN(size(PSTH_ALL1,1)+size(PSTH_ALL2,1),size(PSTH_ALL1,2));
% % % % % % %     PSTH_ALL(1:size(PSTH_ALL1,1),:) = PSTH_ALL1;
% % % % % % %     PSTH_ALL(size(PSTH_ALL1,1)+1:size(PSTH_ALL1,1)+size(PSTH_ALL2,1),:) = PSTH_ALL2;
% % % % % % %     
% % % % % % %     
% % % % % % %     
% % % % % % %     clear PSTH_Z
% % % % % % %     PSTH_Z = (PSTH_ALL);
% % % % % % %     
% % % % % % %     
% % % % % % %     if I==1 PSTH_Z_T = PSTH_Z; end
% % % % % % %     if I==2 PSTH_Z_R = PSTH_Z; end
% % % % % % %     
% % % % % % %     
% % % % % % % end
% % % % % % % 






Signal = spikes.S;
Sigma = 50;

for I=1:2
    
    if I==1 Align_code = 4; end
    if I==2 Align_code = 11; end
    
    if Align_code == 11
        Start = -950; End = 750;
    end
    
    if Align_code == 4
        Start = -450; End = 1250;
    end
    
    
    Align_time = infos(:,Align_code);
   
    
    if I==1 PSTH_T = PSTH_RETURN_n(Signal,Align_time,Start,End,Sigma); end
    if I==2 PSTH_M = PSTH_RETURN_n(Signal,Align_time,Start,End,Sigma); end
    
end


PSTH_Z_T = PSTH_T(:,1+50:size(PSTH_T,2)-50);
PSTH_Z_R = PSTH_M(:,1+50:size(PSTH_M,2)-50);


MIN = nanmin([nanmin(nanmin(PSTH_Z_T)) nanmin(nanmin(PSTH_Z_R))]);
MAX = nanmax([nanmax(nanmax(PSTH_Z_T)) nanmax(nanmax(PSTH_Z_R))]);

clear PSTH_Z_T_temp PSTH_Z_R_temp

PSTH_Z_T_temp(2:size(PSTH_Z_T,1)+1,:) = PSTH_Z_T; 
PSTH_Z_T_temp(1,:)=MIN;
PSTH_Z_T_temp(end+1,:)=MAX;


PSTH_Z_R_temp(2:size(PSTH_Z_R,1)+1,:) = PSTH_Z_R; 
PSTH_Z_R_temp(1,:)=MIN;
PSTH_Z_R_temp(end+1,:)=MAX;


PSTH_Z_T_temp=mat2gray(PSTH_Z_T_temp);
PSTH_Z_R_temp=mat2gray(PSTH_Z_R_temp);

PSTH_Z_T = PSTH_Z_T_temp(2:end-1,:);
PSTH_Z_R = PSTH_Z_R_temp(2:end-1,:);



% smoothing stuff --------------------------------------

span = 20;

TEMP1 = PSTH_Z_T(1:length(Spikes.B),:);
TEMP2 = PSTH_Z_T(length(Spikes.B)+1:length(Spikes.D)+length(Spikes.B),:);

for i=1:size(TEMP1,2)
    TEMP1(:,i) = smooth(TEMP1(:,i),span);
    TEMP2(:,i) = smooth(TEMP2(:,i),span);
end


TEMP3 = PSTH_Z_R(1:length(Spikes.B),:);
TEMP4 = PSTH_Z_R(1+length(Spikes.B):length(Spikes.D)+length(Spikes.B),:);

for i=1:size(TEMP3,2)
    TEMP3(:,i) = smooth(TEMP3(:,i),span);
    TEMP4(:,i) = smooth(TEMP4(:,i),span);
end



PSTH_Z_T_smooth = [TEMP1; TEMP2];
PSTH_Z_R_smooth = [TEMP3; TEMP4];



%% LC and Simple Spike metrics 

FS=8;

clear ff;
ff = figure;

ax1 = subplot(2,2,3);

Start = -400;
End   = 1200;

X = Start:End;
Y = [1:20]';

image(X,Y,flipud(PSTH_Z_T_smooth),'CDataMapping','scaled')
hold on;
plot([0 0],ylim,'-k','Linewidth',1);
plot(xlim, [length(infos)-CHANGE length(infos)-CHANGE],'--k','Linewidth',1);

xlim([Start End]);
set(gca,'yTick',[])
xlabel('Time from target(ms)');
ax1.FontSize=FS;
ylabel([]);
box off;
colormap(jet)
set(gca,'fontsize',8)

% ylim([length(infos)-CHANGE-10 length(infos)-CHANGE+10])
plot([0 0],ylim,'-k','Linewidth',1);
plot(xlim, [length(infos)-CHANGE length(infos)-CHANGE],'--k','Linewidth',1);





ax2 = subplot(2,2,4);
Start = -900;
End   = 700;

X = Start:End;


image(X,Y,flipud(PSTH_Z_R_smooth),'CDataMapping','scaled')
hold on;
plot([0 0],ylim,'-k','Linewidth',1);
plot(xlim, [length(infos)-CHANGE length(infos)-CHANGE],'--k','Linewidth',1);

xlim([Start End]);
set(gca,'yTick',[])
xlabel('Time from movement(ms)');
ax1.FontSize=FS;
ylabel([]);
box off;

colormap(jet)
set(gca,'fontsize',8)

% ylim([length(infos)-CHANGE-20 length(infos)-CHANGE+10])
plot([0 0],ylim,'-k','Linewidth',1);
plot(xlim, [length(infos)-CHANGE length(infos)-CHANGE],'--k','Linewidth',1);



suptitle(NOME)

cd(Results_dir)
filename = strcat(NOME,'_PSTH_SS_JUSTCHANGE');
print(ff, '-dpdf', filename, '-r600')


% print(ff, '-djpeg', 'colorbar', '-r600')





















