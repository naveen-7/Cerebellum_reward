%% Makes the cell plate and gives all necessary information about the neuron/trial behaviour
% Created by NAVEEN ON 06/30/15 at CUMC
% Modified on 04/20/16 at CUMC


% function CER_Cellprops

% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------
% codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
% data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data';


cd(data_dir)
disp('!!! CER_cellprops has started running !!!');

cd(PathName1);
load(filenm);


%% %%%%%%% printing the CELL PLATE %%%%%%%%% %%





% CELL INFORMATION---------------------------------------------


% filename
str = FileName1(6:length(FileName1)-4);
[C,matches] = strsplit(str,{'_'});
FileNAME = str;

% FILE NUMBER --------------------------
FileNumber = C{2};
F_Num = strcat('File # : ', C{2});

% DATE --------------------------
D=C{1}(3:8);
Date_m=D(1:2);
Date_d=D(3:4);
Date_y=strcat('20',D(5:6));
F_Date = strcat('Date : ', Date_y,'-',Date_m,'-',Date_d);
Date = str2num(strcat(Date_y,Date_m,Date_d));

% MONKEY NAME --------------------------
N=C{1}(1:2);
N=upper(N);
if N=='PR'
    F_Monkey='Monkey : Pierre';
elseif N=='BR'
    F_Monkey='Monkey : Barney';
end

Monkey = Monkey_name(1);



%HAND ----------------------------
Hand = WhichHand_n(Spikes.S,CHANGE,Infos);
delete(gcf)
if      Hand==01 F_Hand = 'Hand : Left';
elseif  Hand==02 F_Hand = 'Hand : Right';
elseif  Hand==10 F_Hand = 'Hand : Both';
elseif  Hand==11 F_Hand = 'Hand : Both; L>R';
elseif  Hand==12 F_Hand = 'Hand : Both; R>L';
elseif  Hand==00 F_Hand = 'Hand : Not a hand area';
end




f = figure();
axes('box','off','tickdir','out','Linewidth',1.25,'FontSize',12)

h1=subplot(4,4,2);
axis off;
xl = xlim(h1);
xPos = xl(1) + 0.075;
yl = ylim(h1);
yPos = yl(1) + diff(yl) / 2;

t1 = text(xPos, yPos+0.2,  sprintf(F_Date,   'FontWeight','bold','FontSize',16 ), 'Parent', h1);
t2 = text(xPos, yPos     ,  sprintf(F_Num,    'FontWeight','bold','FontSize',16 ), 'Parent', h1);
t3 = text(xPos, yPos-0.2,  sprintf(F_Monkey, 'FontWeight','bold','FontSize',16 ), 'Parent', h1);
t4 = text(xPos, yPos-0.4 ,  sprintf(F_Hand,   'FontWeight','bold','FontSize',16 ), 'Parent', h1);

CC = get(h1,'position');
% set(h1,'position',[XX(1)-0.12 CC(2) XX(4)/1.5 CC(4)]);


% LEARNING CURVE ---------------------------------------------------
sF = subplot (4,4,3);

hold on;
clear ylim yLim
plot(x_trial_LC,per_corr_LC,'-','color',[1	0.7255	0.0588],'LineWidth',2);
plot(x_trial_LC,per_corr_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
xlabel('Trial number','fontweight','bold','FontSize',7);
ylabel('% Correct','fontweight','bold','FontSize',7);
set(gca,'FontSize',10);


yLim = ylim;
if yLim(1)>20
    ylim([20 yLim(2)]);
end
yLim = ylim;
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);

title('Learning curve','fontweight','bold','FontSize',10);
set(gca,'LineWidth',1)


FF = get(sF,'position'); %l b w h
set(sF,'position',[FF(1) FF(2) FF(3)+0.05 FF(4)]);

% RASTER ----------------------------------------------------
clear xlim;
s3 = subplot (4,4,[6 10]);
title('Simple and Complex Spikes','FontSize',10,'fontweight','bold')
Raster_n(Spikes.S,Infos(:,4),-200,700,[0.6 0.6 0.6],0.15);
xlabel([])
hold on;
Raster_n(Spikes.C,Infos(:,4),-200,700,[0.9333    0.2275    0.5490],2);
xlabel([])
hold on;
plot(xlim,[CHANGE CHANGE],'--k','LineWidth',2);
set(gca,'FontSize',10,'LineWidth',1);
xlabel([])


XX=get(s3,'position'); %l b w h
set(s3,'position',[XX(1)-0.12 XX(2) XX(4)/1.5 XX(4)/1.5]);
set(h1,'position',[XX(1)-0.12 CC(2) XX(4)/1.5 CC(4)]);

% RT -------------------------------------------------------------
% % clear xlim ylim RT_bef RT_aft;
% % n_bin = 0;  % number of bins
% % % w_bin = round(0.1*(size(Infos,1))); % bin width
% % % s_bin = round(0.05*(size(Infos,1)));  % shift in bin
% % 
% % % Part 1: before CHANGE
% % m_bin = round((CHANGE-1)/(s_bin))-1; % max number of bins
% % clear RT_bef 
% % 
% % for i=1:s_bin:CHANGE-1
% %     if(i<=(m_bin-1)*s_bin)
% %         n_bin = n_bin+1;
% %         RT_bef(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
% %     end
% % end
% % 
% % 
% % % Part 2: after CHANGE
% % n_bin = 0;
% % m_bin = round((size(Infos,1)-CHANGE)/(s_bin))-1; % max number of bins
% % clear RT_aft
% % 
% % for i=CHANGE:s_bin:size(Infos,1)
% %     if(i<=((m_bin-2)*s_bin)+CHANGE)
% %         n_bin = n_bin+1;
% %         RT_aft(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
% %     end
% % end
% % 
% % 
% % RT_tot = NaN(length(RT_bef)+length(RT_aft),1);
% % RT_tot(1:length(RT_bef)) = RT_bef;
% % RT_tot(length(RT_bef)+1:length(RT_bef)+length(RT_aft)) = RT_aft;


% RT ------------------------
sA=subplot (4,4,7);
hold on;
plot(x_trial_LC,RT_tot_LC,'-','color',[0.8235    0.4118    0.1176],'LineWidth',2);
plot(x_trial_LC,RT_tot_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
e = errorbar(x_trial_LC,RT_tot_LC,RT_std_tot_LC);
e.Color = [0.8235    0.4118    0.1176];

ylim([min(RT_tot_LC-RT_std_tot_LC) max(RT_tot_LC+RT_std_tot_LC)])
xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
xlabel('Trial number','FontSize',7)
ylabel('RT in ms','FontSize',7)
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);
set(gca,'FontSize',10,'LineWidth',1)
title('RT Profile','FontSize',10,'fontweight','bold')

AA = get(sA,'position'); %l b w h
set(sA,'position',[AA(1) AA(2) AA(3)+0.05 AA(4)]);

% HAND -------------------------------------------------------------
sB = subplot (4,4,11);
Align_time = Infos(5:CHANGE-5,11);
Start_time = -500;
End_time = 600;
Sigma = 15;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.S(5:CHANGE-5);
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
P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1);
P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2);
ylim_L = [min(P_Left) max(P_Left)];
ylim_R = [min(P_Right) max(P_Right)];
% hl = legend('Left','','Right','Location','best');
% set(hl,'FontSize',7);
% legend boxoff\
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
ylabel('Firing Rate Hz','FontSize',7)
xlabel('time in ms','FontSize',7)
if min(ylim_L(1),ylim_R(1))-10 >=10
    ylim([min(ylim_L(1),ylim_R(1))-10 max(ylim_L(2),ylim_R(2))+10]);
else
    ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
end
set(gca,'FontSize',10,'LineWidth',1)
title('Hand Response','fontweight','bold','FontSize',10);

BB = get(sB,'position'); %l b w h
set(sB,'position',[BB(1) BB(2) BB(3)+0.05 BB(4)]);



% MOVEMENT ------------------------------------------------------
sD = subplot (4,4,15);
P = PSTH_n(Spikes.S,Infos(:,11),-500,600,10,[0.2353    0.7020    0.4431]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
title('PSTH @ movement','FontSize',10)
xlabel('time in ms','FontSize',7)
ylabel('Firing Rate Hz','FontSize',7)
if min(P)>=10
    YLIM_2 = [min(P(50:length(P)-50))-5 max(P(50:length(P)-50))+5];
else 
    YLIM_2 = [min(P(50:length(P)-50)) max(P(50:length(P)-50))+2];
end
set(gca,'FontSize',10,'LineWidth',1)
% set(gca,'fontweight','bold')

DD = get(sD,'position'); %l b w h
set(sD,'position',[DD(1) DD(2) DD(3)+0.05 DD(4)]);

MOV_max = nanmax(P);
MOV_min = nanmin(P);

% repositioning the raster -----------

% subplot(4,4,[6 10])
Ynew = 0.05 + [ (BB(4)*2) + (AA(2)-BB(2)-BB(4)) - (XX(4)) ] / 2;
set(s3,'position',[XX(1)-0.12 Ynew+BB(2) XX(4)/1.5 XX(4)/1.5]);  %l b w h


% TARGET ------------------------------------------------------
s7 = subplot (4,4,14);
P = PSTH_n(Spikes.S,Infos(:,4),-200,700,10,[0.2118    0.3922    0.5451]);
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
title(' PSTH @ target','FontSize',10)
xlabel('time in ms','FontSize',7)
ylabel('Firing Rate Hz','FontSize',7)
if min(P)>=10
    YLIM_1 = [min(P(50:length(P)-50))-5 max(P(50:length(P)-50))+5];
else 
    YLIM_1 = [min(P(50:length(P)-50)) max(P(50:length(P)-50))+2];
end
YY = get(s7,'position');
set(s7,'position',[XX(1)-0.12 YY(2) XX(4)/1.5 YY(4)]);
set(gca,'FontSize',10,'LineWidth',1)
% set(gca,'fontweight','bold')
hold on;



VIS_max = nanmax(P);
VIS_min = nanmin(P);

linkaxes([s7,s3],'x');


set(s7,'ylim',([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]));
% hold on;
% plot([0 0],ylim,'-k');
% subplot(4,4,15)
set(sD,'ylim',([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]));


filename = strcat(FileName1(6:17),'_','props');
% set(gcf, 'PaperUnits','inches','PaperSize',[8.5 11],'PaperPosition',[-2 -0.5 12 9])
set(gcf, 'PaperUnits','inches','PaperSize',[8.5 11],'PaperPosition',[-1.2 0 12 9])
cd(Results_dir);
print(f, '-dpdf', filename, '-r400')
%delete(gcf)


filename_cell = filename;








MAX_RATE = nanmax([VIS_max MOV_max])
MIN_RATE = nanmin([VIS_min MOV_min])

cd(PoP_dir);
save(POP_file,'MAX_RATE','MIN_RATE','-append');
save(DATA_file,'MAX_RATE','MIN_RATE','-append');











% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% Another cell plate but only main task
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % F= figure();
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % CELL INFORMATION---------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % clear CC matches;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % filename
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % str = FileName1(6:length(FileName1)-4);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % [C,matches] = strsplit(str,{'_'});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % FileNAME = str;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % FILE NUMBER --------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % FileNumber = C{2};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % F_Num = strcat('File # : ', C{2});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hh = suptitle(strcat(C{1},'-',C{2}));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(hh,'FontSize',10,'FontWeight','bold')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % LEARNING CURVE ---------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot (3,2,1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % clear ylim yLim
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot(x_trial_LC,per_corr_LC,'-','color',[1	0.7255	0.0588],'LineWidth',2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot(x_trial_LC,per_corr_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % xlabel('Trial number','fontweight','bold','FontSize',7);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % ylabel('% Correct','fontweight','bold','FontSize',7);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(gca,'FontSize',10);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % yLim = ylim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % if yLim(1)>20
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim([20 yLim(2)]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % yLim = ylim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % title('Learning curve','fontweight','bold','FontSize',10);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(gca,'LineWidth',1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % RT ------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot (3,2,2)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot(x_trial_LC,RT_tot_LC,'-','color',[0.8235    0.4118    0.1176],'LineWidth',2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot(x_trial_LC,RT_tot_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % ylim([min(RT_tot_LC-10) max(RT_tot_LC+10)])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % xlabel('Trial number','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % ylabel('RT in ms','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(gca,'FontSize',10,'LineWidth',1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % title('RT Profile','FontSize',10,'fontweight','bold')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % RASTER @T----------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     clear xlim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,3);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         title('@ target','FontSize',8,'fontweight','bold');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Raster_n(Spikes.S(1:CHANGE),Infos(1:CHANGE,4),-600,1200,[0.6 0.6 0.6],0.5,0,1,'o');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Raster_n(Spikes.C(1:CHANGE),Infos(1:CHANGE,4),-600,1200,[0.9333    0.2275    0.5490],1.5,0,1,'o');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'FontSize',7,'LineWidth',3);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontweight','bold')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlabel([])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     axis off;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontsize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % RASTER @M----------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     clear xlim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,4);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         title('@ movement','FontSize',8,'fontweight','bold');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Raster_n(Spikes.S(1:CHANGE),Infos(1:CHANGE,11),-800,700,[0.6 0.6 0.6],0.5,0,1,'o');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Raster_n(Spikes.C(1:CHANGE),Infos(1:CHANGE,11),-800,700,[0.9333    0.2275    0.5490],1.5,0,1,'o');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'FontSize',7,'LineWidth',3);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontweight','bold')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlabel([])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     axis off;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontsize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % HAND_@T -------------------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,5);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Align_time = Infos(1:CHANGE,4);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Start_time = -650;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     End_time = 1250;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sigma = 20;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Colour1 = [0.2 0.2 0.2]; % DARK,  Left
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Colour2 = [0.7 0.7 0.7]; % LIGHT, Right
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     LW = 1.25;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     clear SIGNAL xlim xLIM;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     SIGNAL = Spikes.S(1:CHANGE);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sig_L = SIGNAL;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sig_R = SIGNAL;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     for i=1:length(SIGNAL)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if (Infos(i,9)==1)      % Left
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Sig_L{i,1}=[];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if (Infos(i,9)==0)      % Right
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Sig_R{i,1}=[];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim_L = [min(P_Left(50:length(P_Left)-50)) max(P_Left(50:length(P_Left)-50))];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim_R = [min(P_Right(50:length(P_Right)-50)) max(P_Right(50:length(P_Right)-50))];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xLIM=xlim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlim([xLIM(1)+50 xLIM(2)-50]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylabel('Firing Rate Hz','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlabel('time in ms','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if min(ylim_L(1),ylim_R(1))-10 >=10
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ylim([min(ylim_L(1),ylim_R(1))-10 max(ylim_L(2),ylim_R(2))+10]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         else
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'FontSize',7,'LineWidth',1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylabel('Sp/s');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontsize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     YLIM_T = ylim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     clear ylim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % HAND_@R -------------------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,6);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Align_time = Infos(1:CHANGE,11);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Start_time = -850;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     End_time = 750;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sigma = 20;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Colour1 = [0.2 0.2 0.2];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Colour2 = [0.7 0.7 0.7];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     clear SIGNAL xlim xLIM;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     SIGNAL = Spikes.S(1:CHANGE);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sig_L = SIGNAL;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Sig_R = SIGNAL;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     for i=1:length(SIGNAL)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if (Infos(i,9)==1)      % Left
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Sig_L{i,1}=[];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if (Infos(i,9)==0)      % Right
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Sig_R{i,1}=[];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     P_Left  = PSTH_n(Sig_L,Align_time,Start_time,End_time,Sigma,Colour1,LW);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     P_Right = PSTH_n(Sig_R,Align_time,Start_time,End_time,Sigma,Colour2,LW);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim_L = [min(P_Left(50:length(P_Left)-50)) max(P_Left(50:length(P_Left)-50))];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim_R = [min(P_Right(50:length(P_Right)-50)) max(P_Right(50:length(P_Right)-50))];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xLIM=xlim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlim([xLIM(1)+50 xLIM(2)-50]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylabel('Firing Rate Hz','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlabel('time in ms','FontSize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         if min(ylim_L(1),ylim_R(1))-10 >=10
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ylim([min(ylim_L(1),ylim_R(1))-10 max(ylim_L(2),ylim_R(2))+10]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         else
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'FontSize',7,'LineWidth',1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylabel('Sp/s');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'fontsize',7)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     YLIM_M = ylim;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,5)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     subplot (3,2,6)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     ylim([min(YLIM_T(1),YLIM_M(1)) max(YLIM_T(2),YLIM_M(2))]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hold on;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     disp(' >>> Almost DONE <<<');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % filename = strcat(FileName1(6:17),'_','props_naturale');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(gcf, 'PaperUnits','inches','PaperSize',[8.5 11],'PaperPosition',[-1.2 0 12 9])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % cd(Results_dir);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % print(F, '-dpdf', filename, '-r400')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %delete(gcf)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% %% Another cell plate but for all trials

F= figure();
% CELL INFORMATION---------------------------------------------

clear CC matches;
% filename
str = FileName1(6:length(FileName1)-4);
[C,matches] = strsplit(str,{'_'});
FileNAME = str;

% FILE NUMBER --------------------------
FileNumber = C{2};
F_Num = strcat('File # : ', C{2});


hh = suptitle(strcat(C{1},'-',C{2}));
set(hh,'FontSize',10,'FontWeight','bold')


% LEARNING CURVE ---------------------------------------------------



subplot (3,2,1)

hold on;
clear ylim yLim
plot(x_trial_LC,per_corr_LC,'-','color',[1	0.7255	0.0588],'LineWidth',2);
plot(x_trial_LC,per_corr_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
xlabel('Trial number','fontweight','bold','FontSize',7);
ylabel('% Correct','fontweight','bold','FontSize',7);
set(gca,'FontSize',10);


yLim = ylim;
if yLim(1)>20
    ylim([20 yLim(2)]);
end
yLim = ylim;
plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);

title('Learning curve','fontweight','bold','FontSize',10);
set(gca,'LineWidth',1)


%
% % RT ------------------------
% subplot (3,2,2)
% hold on;
% plot(x_trial_LC,RT_tot_LC,'-','color',[0.8235    0.4118    0.1176],'LineWidth',2);
% plot(x_trial_LC,RT_tot_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
%
% ylim([min(RT_tot_LC-10) max(RT_tot_LC+10)])
% xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
% xlabel('Trial number','FontSize',7)
% ylabel('RT in ms','FontSize',7)
% plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
% plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);
% set(gca,'FontSize',10,'LineWidth',1)
% title('RT Profile','FontSize',10,'fontweight','bold')




% Start_time_T = -650;
% End_time_T = 1750;
% Start_time_R = -1050;
% End_time_R = 750;


Start_time_T = -450;
End_time_T = 1250;
Start_time_R = -950;
End_time_R = 750;





% HAND_@R -------------------------------------------------------------
subplot (3,2,2);
Align_time = Infos(1:CHANGE,11);
Start_time = -850;
End_time = 750;
Sigma = 20;
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
clear SIGNAL xlim xLIM;
SIGNAL = Spikes.S(1:CHANGE);
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
LW=1;
P_Left  = PSTH_n(Sig_L,Align_time,Start_time_R,End_time_R,Sigma,Colour1,LW);
P_Right = PSTH_n(Sig_R,Align_time,Start_time_R,End_time_R,Sigma,Colour2,LW);
ylim_L = [min(P_Left) max(P_Left)];
ylim_R = [min(P_Right) max(P_Right)];
xLIM=xlim;
xlim([xLIM(1)+50 xLIM(2)-50]);
ylabel('Firing Rate Hz','FontSize',7)
xlabel('time in ms','FontSize',7)
if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
    if min(ylim_L(1),ylim_R(1))-10 >=10
        ylim([min(ylim_L(1),ylim_R(1))-10 max(ylim_L(2),ylim_R(2))+10]);
    else
        ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
    end
end
set(gca,'FontSize',7,'LineWidth',1)
ylabel([]);
set(gca,'fontsize',7)

YLIM_M = ylim;




% RASTER @T----------------------------------------------------
clear xlim;
subplot (3,2,3);
title('@ target','FontSize',8,'fontweight','bold');
Raster_n(Spikes.S,Infos(:,4),Start_time_T+50,End_time_T-50,[0.6 0.6 0.6],0.25,0,1,'o');
hold on;
Raster_n(Spikes.C,Infos(:,4),Start_time_T+50,End_time_T-50,[0.9333    0.2275    0.5490],1,0,1,'o');
set(gca,'FontSize',7,'LineWidth',3);
set(gca,'fontweight','bold')
xlabel([])
axis off;
set(gca,'fontsize',7)
hold on;
plot(xlim,[CHANGE CHANGE],'-k','linewidth',1);
plot(xlim,[LEARNT LEARNT],'-k','linewidth',1);


% RASTER @M----------------------------------------------------
clear xlim;
subplot (3,2,4);
title('@ movement','FontSize',8,'fontweight','bold');
Raster_n(Spikes.S,Infos(:,11),Start_time_R+50,End_time_R-50,[0.6 0.6 0.6],0.25,0,1,'o');
hold on;
Raster_n(Spikes.C,Infos(:,11),Start_time_R+50,End_time_R-50,[0.9333    0.2275    0.5490],1,0,1,'o');
set(gca,'FontSize',7,'LineWidth',3);
set(gca,'fontweight','bold')
xlabel([])
axis off;
set(gca,'fontsize',7)
hold on;
plot(xlim,[CHANGE CHANGE],'-k','linewidth',1);
plot(xlim,[LEARNT LEARNT],'-k','linewidth',1);


% TARGET ------------------------------------------------------

% limit = round(0.25*(LEARNT-CHANGE-5));

Bef_Sig = Spikes.S(CHANGE-9:CHANGE);
Dur_Sig = Spikes.S(CHANGE+5:CHANGE+15);
Aft_Sig = Spikes.S(length(Infos)-9:length(Infos));



subplot (3,2,5);
Sigma = 10;
hold on;
cd(codes_dir);
P1 = PSTH_n(Bef_Sig,Infos((CHANGE-9:CHANGE),4),Start_time_T,End_time_T,20,[0.2118    0.3922    0.5451]); % Before
P2 = PSTH_n(Dur_Sig,Infos((CHANGE+5:CHANGE+15),4),Start_time_T,End_time_T,20,[0.2353    0.7020    0.4431]); % During
P3 = PSTH_n(Aft_Sig,Infos((length(Infos)-9:length(Infos)),4),Start_time_T,End_time_T,20,[0.9333    0.6039         0]); % After
xLim = xlim;
xlim([xLim(1)+50 xLim(2)-50]);
set(gca,'linewidth',1)
xlabel('time in ms','FontSize',7)
ylabel('Sp/s','FontSize',7)
set(gca,'fontsize',7)
MAXMAX = nanmax([nanmax(P1) nanmax(P2) nanmax(P3)]);
MINMIN = nanmin([nanmin(P1(50:length(P1)-50)) nanmin(P2(50:length(P2)-50)) nanmin(P3(50:length(P3)-50))]);

if MINMIN>=10
    YLIM_1 = [MINMIN-5 MAXMAX+5];
else
    YLIM_1 = [MINMIN MAXMAX+2];
end


% MOVEMENT ------------------------------------------------------

subplot(3,2,6)
hold on;

P1 = PSTH_n(Bef_Sig,Infos((CHANGE-9:CHANGE),11),Start_time_R,End_time_R,20,[0.2118    0.3922    0.5451]); % Before
P2 = PSTH_n(Dur_Sig,Infos((CHANGE+5:CHANGE+15),11),Start_time_R,End_time_R,20,[0.2353    0.7020    0.4431]); % During
P3 = PSTH_n(Aft_Sig,Infos((length(Infos)-9:length(Infos)),11),Start_time_R,End_time_R,20,[0.9333    0.6039         0]); % After
xLim = xlim;
xlim([xLim(1)+50 xLim(2)-50]);
set(gca,'linewidth',1)
xlabel('time in ms','FontSize',7)
ylabel('Sp/s','FontSize',7)
set(gca,'fontsize',7)
MAXMAX = nanmax([nanmax(P1) nanmax(P2) nanmax(P3)]);
MINMIN = nanmin([nanmin(P1(50:length(P1)-50)) nanmin(P2(50:length(P2)-50)) nanmin(P3(50:length(P3)-50))]);

if MINMIN>=10
    YLIM_2 = [MINMIN-5 MAXMAX+5];
else
    YLIM_2 = [MINMIN MAXMAX+2];
end



subplot (3,2,5)
ylim([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
subplot (3,2,6)
ylim([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);




disp(' >>> Almost DONE <<<');


filename = strcat(FileName1(6:17),'_','props_ALL');
% set(gcf, 'PaperUnits','inches','PaperSize',[8.5 11],'PaperPosition',[-1.2 0 12 9])
cd(Results_dir);
print(F, '-dpdf', filename, '-r400')
%delete(gcf)










disp('!!! END OF CODE CER_cellprops !!!');

% end