




cd(data_dir)
disp('!!! Code has started running !!!');

cd(PathName1);
load(filenm);


for i=2:length(Infos)
    
    if (Infos(i-1,10)==2)  
        OUTCOME(i,1)=2;
    end
    
    if (Infos(i-1,10)==1)     % Correct 
        OUTCOME(i,1)=1;
    end
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% % 
% % 
% % subplot (3,2,1)
% % 
% % hold on;
% % clear ylim yLim
% % plot(x_trial_LC,per_corr_LC,'-','color',[1	0.7255	0.0588],'LineWidth',2);
% % plot(x_trial_LC,per_corr_LC,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% % xlabel('Trial number','fontweight','bold','FontSize',7);
% % ylabel('% Correct','fontweight','bold','FontSize',7);
% % set(gca,'FontSize',10);
% % 
% % 
% % yLim = ylim;
% % if yLim(1)>20
% %     ylim([20 yLim(2)]);
% % end
% % yLim = ylim;
% % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color',[0.4 0.4 0.4]);
% % xlim([x_trial_LC(1) x_trial_LC(length(x_trial_LC))]);
% % plot([LEARNT LEARNT], ylim,'--','lineWidth',1,'color',[0.7 0.7 0.7]);
% % 
% % title('Learning curve','fontweight','bold','FontSize',10);
% % set(gca,'LineWidth',1)
% % 
% % 
% % Start_time_T = -450;
% % End_time_T = 1250;
% % Start_time_R = -950;
% % End_time_R = 750;
% % 
% % 
% % 
% % 
% % 
% % HAND_@R -------------------------------------------------------------
% % subplot (3,2,2);
% % Align_time = Infos(1:CHANGE,11);
% % Start_time = -850;
% % End_time = 750;
% % Sigma = 20;
% % Colour1 = [0.2 0.2 0.2];
% % Colour2 = [0.7 0.7 0.7];
% % clear SIGNAL xlim xLIM;
% % SIGNAL = Spikes.S(1:CHANGE);
% % Sig_L = SIGNAL;
% % Sig_R = SIGNAL;
% % 
% % for i=1:length(SIGNAL)
% %     if (Infos(i,9)==1)      % Left
% %         Sig_L{i,1}=[];
% %     end
% %     if (Infos(i,9)==0)      % Right
% %         Sig_R{i,1}=[];
% %     end
% % end
% % LW=1;
% % P_Left  = PSTH_n(Sig_L,Align_time,Start_time_R,End_time_R,Sigma,Colour1,LW);
% % P_Right = PSTH_n(Sig_R,Align_time,Start_time_R,End_time_R,Sigma,Colour2,LW);
% % ylim_L = [min(P_Left) max(P_Left)];
% % ylim_R = [min(P_Right) max(P_Right)];
% % xLIM=xlim;
% % xlim([xLIM(1)+50 xLIM(2)-50]);
% % ylabel('Firing Rate Hz','FontSize',7)
% % xlabel('time in ms','FontSize',7)
% % if ~isnan(ylim_L(1)) & ~isnan(ylim_L(2)) & ~isnan(ylim_R(1)) & ~isnan(ylim_R(2))
% %     if min(ylim_L(1),ylim_R(1))-10 >=10
% %         ylim([min(ylim_L(1),ylim_R(1))-10 max(ylim_L(2),ylim_R(2))+10]);
% %     else
% %         ylim([min(ylim_L(1),ylim_R(1)) max(ylim_L(2),ylim_R(2))+2]);
% %     end
% % end
% % set(gca,'FontSize',7,'LineWidth',1)
% % ylabel([]);
% % set(gca,'fontsize',7)
% % 
% % YLIM_M = ylim;
% % 
% % 
% % 
% % 
% % RASTER @T----------------------------------------------------
% % clear xlim;
% % subplot (3,2,3);
% % title('@ target','FontSize',8,'fontweight','bold');
% % Raster_n(Spikes.S,Infos(:,4),Start_time_T+50,End_time_T-50,[0.6 0.6 0.6],0.25,0,1,'o');
% % hold on;
% % Raster_n(Spikes.C,Infos(:,4),Start_time_T+50,End_time_T-50,[0.9333    0.2275    0.5490],1,0,1,'o');
% % set(gca,'FontSize',7,'LineWidth',3);
% % set(gca,'fontweight','bold')
% % xlabel([])
% % axis off;
% % set(gca,'fontsize',7)
% % hold on;
% % plot(xlim,[CHANGE CHANGE],'-k','linewidth',1);
% % plot(xlim,[LEARNT LEARNT],'-k','linewidth',1);
% % 
% % 
% % RASTER @M----------------------------------------------------
% % clear xlim;
% % subplot (3,2,4);
% % title('@ movement','FontSize',8,'fontweight','bold');
% % Raster_n(Spikes.S,Infos(:,11),Start_time_R+50,End_time_R-50,[0.6 0.6 0.6],0.25,0,1,'o');
% % hold on;
% % Raster_n(Spikes.C,Infos(:,11),Start_time_R+50,End_time_R-50,[0.9333    0.2275    0.5490],1,0,1,'o');
% % set(gca,'FontSize',7,'LineWidth',3);
% % set(gca,'fontweight','bold')
% % xlabel([])
% % axis off;
% % set(gca,'fontsize',7)
% % hold on;
% % plot(xlim,[CHANGE CHANGE],'-k','linewidth',1);
% % plot(xlim,[LEARNT LEARNT],'-k','linewidth',1);
% % 

% TARGET ------------------------------------------------------

% limit = round(0.25*(LEARNT-CHANGE-5));

Start_time_T = -450;
End_time_T = 1250;
Start_time_R = -950;
End_time_R = 750;

for ii=1:2
    
    IND = find(OUTCOME==ii);
BEF_IND = find(OUTCOME(CHANGE-9:CHANGE)==ii)+CHANGE-9-1;
DUR_IND = find(OUTCOME(CHANGE+5:CHANGE+15)==ii)+CHANGE+5-1;
AFT_IND = find(OUTCOME(length(Infos)-9:length(Infos))==ii)+length(Infos)-9-1;

Bef_Sig = Spikes.S(BEF_IND);
Dur_Sig = Spikes.S(DUR_IND);
Aft_Sig = Spikes.S(AFT_IND);



subplot (3,2,2*(ii-1)+1+2);
Sigma = 10;
hold on;
cd(codes_dir);
P1 = PSTH_n(Bef_Sig,Infos((BEF_IND),4),Start_time_T,End_time_T,20,[0.2118    0.3922    0.5451]); % Before
P2 = PSTH_n(Dur_Sig,Infos((DUR_IND),4),Start_time_T,End_time_T,20,[0.2353    0.7020    0.4431]); % During
P3 = PSTH_n(Aft_Sig,Infos((AFT_IND),4),Start_time_T,End_time_T,20,[0.9333    0.6039         0]); % After
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

subplot(3,2,2*(ii-1)+1+3)
hold on;

P1 = PSTH_n(Bef_Sig,Infos((BEF_IND),11),Start_time_R,End_time_R,20,[0.2118    0.3922    0.5451]); % Before
P2 = PSTH_n(Dur_Sig,Infos((DUR_IND),11),Start_time_R,End_time_R,20,[0.2353    0.7020    0.4431]); % During
P3 = PSTH_n(Aft_Sig,Infos((AFT_IND),11),Start_time_R,End_time_R,20,[0.9333    0.6039         0]); % After
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



subplot (3,2,2*(ii-1)+1+2)
ylim([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);
subplot (3,2,2*(ii-1)+1+3)
ylim([min(YLIM_1(1),YLIM_2(1)) max(YLIM_1(2),YLIM_2(2))]);
hold on;
plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',1);

clear YLIM_1 YLIM_2

end

disp(' >>> Almost DONE <<<');

subplot (3,2,1)
axis off;
text(0,0.2,'top: prev correct')
text(0,0,'bottom: prev wrong')

subplot (3,2,2)
axis off;
text(0,0.4,'Before','color',[0.2118    0.3922    0.5451])
text(0,0.2,'During','color',[0.2353    0.7020    0.4431])
text(0,0,'After','color',[0.9333    0.6039         0])

filename = strcat(FileName1(6:17),'_','props_CORR_WRONG');
% set(gcf, 'PaperUnits','inches','PaperSize',[8.5 11],'PaperPosition',[-1.2 0 12 9])
cd(Results_dir);
print(F, '-dpdf', filename, '-r400')
%delete(gcf)




