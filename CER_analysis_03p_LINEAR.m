

%% function CER_analysis_03c_LINEAR
% Memory analysis

NUM = 24;

%%
% 
WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;



for i=2:length(Infos)
    
    if (Infos(i,10)==1) && (Infos(i-1,10)==2)     % Wrong -> Correct
        Infos(i,16)=1;
    end
    
    if (Infos(i,10)==2) && (Infos(i-1,10)==1)     % Correct -> Wrong
        Infos(i,16)=2;
    end
    
    if (Infos(i,10)==1) && (Infos(i-1,10)==1)     % Correct -> Correcct
        Infos(i,16)=3;
    end
    
    if (Infos(i,10)==2) && (Infos(i-1,10)==2)     % Wrong -> Wrong
        Infos(i,16)=4;
    end
    
end



% Start_T = -450;
% End_T  =1050;
% StartLIM_T = -400; EndLIM_T = 800;
% 
% Start_M = -750;
% End_M  =750;
% StartLIM_M = -500; EndLIM_M = 700;


%%% NEW  %%% 11/27/19
Start_T = -1550;
End_T  =1550;
StartLIM_T = -1400; EndLIM_T = 400;

Start_M = -1550;
End_M  =1550;
StartLIM_M = -400; EndLIM_M = 1400;

FF = figure;
SIG =30;


% Wrong -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==1)+CHANGE-1;
WC_NUM=length(IND);
subplot(4,4,1)
WC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
subplot(4,4,2)
WC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim; 


% Correct -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==2)+CHANGE-1;
CW_NUM=length(IND);
subplot(4,4,1)
CW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim;
subplot(4,4,2)
CW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Correct -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==3)+CHANGE-1;
CC_NUM=length(IND);
subplot(4,4,1)
CC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
subplot(4,4,2)
CC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;

% Wrong -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==4)+CHANGE-1;
WW_NUM=length(IND);
subplot(4,4,1)
WW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim;
subplot(4,4,2)
WW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;


%MAIN TASK
subplot(4,4,1)
IND = CHANGE-15:CHANGE;
MT_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim;
subplot(4,4,2)
MT_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM10 = ylim; 


subplot(4,4,1)
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
hold on;
plot([0 0],ylim,'-k','linewidth',1);
subplot(4,4,2) 
hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
plot([0 0],ylim,'-k','linewidth',1); 



clear temp1 temp2 temp3 temp4; 
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==4)+CHANGE-1;
temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); 
temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); 
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==2)+CHANGE-1;
temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG);
temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG);
WN_T = [temp1;temp2]; WN_M = [temp3;temp4]; WN = [WN_T WN_M];

clear temp1 temp2 temp3 temp4; 
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==1)+CHANGE-1;
temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); 
temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); 
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==3)+CHANGE-1;
temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG);
temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG);
CN_T = [temp1;temp2]; CN_M = [temp3;temp4]; CN = [CN_T CN_M];

clear P_ttest
parfor ii=1:size(WN,2)
    P_ttest(1,ii) = ttest_NN(WN(:,ii),CN(:,ii));
end



subplot(4,4,9); hold on;
errorline_n(Start_T:End_T,nanmean(WN_T),nanstd(WN_T)/sqrt(size(WN_T,1)),1,[1 0 0]);
errorline_n(Start_T:End_T,nanmean(CN_T),nanstd(CN_T)/sqrt(size(CN_T,1)),1,[0 0 1]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IND1 = CHANGE-20:CHANGE;
IND2 = CHANGE:CHANGE+20;
IND3 = CHANGE+21:CHANGE+40;
IND4 = CHANGE+40:nanmin([LEARNT CHANGE+60 size(Infos,1)]);

IND1 = IND1(find(Infos(IND1,10)==1));
IND2 = IND2(find(Infos(IND2,10)==1));
IND3 = IND3(find(Infos(IND3,10)==1));
IND4 = IND4(find(Infos(IND4,10)==1));


subplot(4,4,3); hold on;
PSTH_n(Spikes.S(IND1,:),Infos(IND1,4),Start_T,End_T,SIG,'k',1,0); % Aligned to target
PSTH_n(Spikes.S(IND2,:),Infos(IND2,4),Start_T,End_T,SIG,[0.3 0.3 0.3],1,0); % Aligned to target
PSTH_n(Spikes.S(IND3,:),Infos(IND3,4),Start_T,End_T,SIG,[0.6 0.6 0.6],1,0); % Aligned to target
PSTH_n(Spikes.S(IND4,:),Infos(IND4,4),Start_T,End_T,SIG,[0.9 0.9 0.9],1,0); hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
plot([0 0],ylim,'-k','linewidth',1); xlim([StartLIM_T EndLIM_T]);


subplot(4,4,4); hold on;
PSTH_n(Spikes.S(IND1,:),Infos(IND1,11),Start_M,End_M,SIG,'k',1,0); % Aligned to target
PSTH_n(Spikes.S(IND2,:),Infos(IND2,11),Start_M,End_M,SIG,[0.3 0.3 0.3],1,0); % Aligned to target
PSTH_n(Spikes.S(IND3,:),Infos(IND3,11),Start_M,End_M,SIG,[0.6 0.6 0.6],1,0); % Aligned to target
PSTH_n(Spikes.S(IND4,:),Infos(IND4,11),Start_M,End_M,SIG,[0.9 0.9 0.9],1,0); hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
plot([0 0],ylim,'-k','linewidth',1); xlim([StartLIM_M EndLIM_M]); set(gca,'YTick',[]); ylabel([]); xlabel([]);


P1_T = PSTH_RETURN_n(Spikes.S(IND1,:),Infos(IND1,4),Start_T,End_T,SIG); % Aligned to target
P2_T = PSTH_RETURN_n(Spikes.S(IND2,:),Infos(IND2,4),Start_T,End_T,SIG); % Aligned to target
P3_T = PSTH_RETURN_n(Spikes.S(IND3,:),Infos(IND3,4),Start_T,End_T,SIG); % Aligned to target
P4_T = PSTH_RETURN_n(Spikes.S(IND4,:),Infos(IND4,4),Start_T,End_T,SIG); % Aligned to target

P1_M = PSTH_RETURN_n(Spikes.S(IND1,:),Infos(IND1,11),Start_M,End_M,SIG); % Aligned to target
P2_M = PSTH_RETURN_n(Spikes.S(IND2,:),Infos(IND2,11),Start_M,End_M,SIG); % Aligned to target
P3_M = PSTH_RETURN_n(Spikes.S(IND3,:),Infos(IND3,11),Start_M,End_M,SIG); % Aligned to target
P4_M = PSTH_RETURN_n(Spikes.S(IND4,:),Infos(IND4,11),Start_M,End_M,SIG); % Aligned to target


P1 = [P1_T P1_M]; P2 = [P2_T P2_M]; P3 = [P3_T P3_M]; P4 = [P4_T P4_M];
LENN = nanmax([size(P1,1) size(P2,1) size(P3,1) size(P4,1)]);


parfor ii=1:size(P1,2)
   ANOVA_MAT = NaN(LENN,4);
   ANOVA_MAT(1:length(P1(:,ii)),1) = P1(:,ii);
   ANOVA_MAT(1:length(P2(:,ii)),2) = P2(:,ii);
   ANOVA_MAT(1:length(P3(:,ii)),3) = P3(:,ii);
   ANOVA_MAT(1:length(P4(:,ii)),4) = P4(:,ii);
   P_anova(1,ii) = anova1(ANOVA_MAT);
   delete(gcf); delete(gcf);
end

alltime = [Start_T:End_T Start_M:End_M];
LIM1 = find(alltime==StartLIM_T,1);
LIM2 = find(alltime==EndLIM_T,1);
LIM3 = find(alltime==StartLIM_M); LIM3=LIM3(2);
LIM4 = find(alltime==EndLIM_M); LIM4 = LIM4(2);

S1 = subplot(4,4,[5]); hold on;
plot(StartLIM_T:EndLIM_T,log10(P_ttest(LIM1:LIM2)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_T EndLIM_T]);
S2 = subplot(4,4,[6]); hold on;
plot(StartLIM_M:EndLIM_M,log10(P_ttest(LIM3:LIM4)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_M EndLIM_M]);


S3 = subplot(4,4,[7]); hold on;
plot(StartLIM_T:EndLIM_T,log10(P_anova(LIM1:LIM2)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_T EndLIM_T]);
S4 = subplot(4,4,[8]); hold on;
plot(StartLIM_M:EndLIM_M,log10(P_anova(LIM3:LIM4)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_M EndLIM_M]);


linkaxes([S1 S2 S3 S4],'y');


P_anova(isnan(P_anova))=1;
P_ttest(isnan(P_ttest))=1;
CORR_VAL = corr(P_anova(:),P_ttest(:));

STATE_DELTA.P_anova = P_anova;
STATE_DELTA.P_ttest = P_ttest;
STATE_DELTA.CORR = CORR_VAL;

suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = 'PSTH_Simple spike_CUATRO_STATECHANGE';
print(FF, '-dpdf', filename, '-r400')


save(MERGE_file,'STATE_DELTA','-append');
save(ALLCELLS_file,'STATE_DELTA','-append');
save(POP_file,'STATE_DELTA','-append');


% end
