
%% function ZERO_manipulandam_NUEVO2

% Written by naveen on 1/5/18 at cumc

% Take 5 trials before change and 5 trials after change
% take 3 random trials in bef and 3 random trials in aft and take mean
% repeat 500 times
% Randomize the trials and repeat

Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\MERGED\',NOME_COMBINED);
cd(Results_dir)


% NUM=round(52/2);
NUM=50;


HAND_ONE = [ONE.H; ONE.V];
HAND_TWO = [TWO.H; TWO.V];

HAND_ONE_H = ONE.H;
HAND_ONE_V = ONE.V;
HAND_TWO_H = TWO.H;
HAND_TWO_V = TWO.V;








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % for inCOUNT=1:50



% for inCOUNT=1:20


clear real_H real_V;

% ACROSS GROUPS ------------------- H
for i=1:NUM
        clear TEMP1 TEMP2;
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_H(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_H(:,ind2),2);
    
%     real_H(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    real_H(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
end


TEMP1_across_H = TEMP1;
TEMP2_across_H = TEMP2;



% ACROSS GROUPS ------------------- V
for i=1:NUM
        clear TEMP1 TEMP2;
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_V(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_V(:,ind2),2);
    
%     real_V(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    real_V(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
end


% real = nanmean([real_H; real_V],1);
% real(find(real>nanmean(real)+2*nanstd(real)))=NaN;



TEMP1_across_V = TEMP1;
TEMP2_across_V = TEMP2;








clear shuffle_H shuffle_V;


% WITHIN GROUPS ------------------- H
for i=1:NUM/2
        clear TEMP1 TEMP2;
    ind = randperm(5); 
    ind1=ind(1:2);
    ind2=ind(3:4);
    TEMP1 = nanmean(HAND_ONE_H(:,ind1),2);
    TEMP2 = nanmean(HAND_ONE_H(:,ind2),2);
    
%     shuffle_H(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    shuffle_H(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
end


TEMP1_H1 = TEMP1;
TEMP2_H1 = TEMP2;

TEMP1_within1_H1 = TEMP1;
TEMP2_within1_H1 = TEMP2;


% WITHIN GROUPS ------------------- V
for i=1:NUM/2
        clear TEMP1 TEMP2;
    ind = randperm(5); 
    ind1=ind(1:2);
    ind2=ind(3:4);
    TEMP1 = nanmean(HAND_ONE_V(:,ind1),2);
    TEMP2 = nanmean(HAND_ONE_V(:,ind2),2);
    
%     shuffle_V(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    shuffle_V(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);

end


TEMP1_V1 = TEMP1;
TEMP2_V1 = TEMP2;

TEMP1_within1_V = TEMP1;
TEMP2_within1_V = TEMP2;


%%

% WITHIN GROUPS ------------------- H
for i=NUM/2+1:NUM
        clear TEMP1 TEMP2;
    ind = randperm(5); 
    ind1=ind(1:2);
    ind2=ind(3:4);
    TEMP1 = nanmean(HAND_TWO_H(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_H(:,ind2),2);
    
%     shuffle_H(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    shuffle_H(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
end


TEMP1_H2 = TEMP1;
TEMP2_H2 = TEMP2;

TEMP1_within1_H2 = TEMP1;
TEMP2_within1_H2 = TEMP2;



TEMP1_H = nanmean([TEMP1_H1, TEMP1_H2],2);
TEMP2_H = nanmean([TEMP2_H1, TEMP2_H2],2);

TEMP1_within1_H = nanmean([TEMP1_H1, TEMP1_H2],2);
TEMP2_within1_H = nanmean([TEMP2_H1, TEMP2_H2],2);





% WITHIN GROUPS ------------------- V
for i=NUM/2+1:NUM
        clear TEMP1 TEMP2;
    ind = randperm(5); 
    ind1=ind(1:2);
    ind2=ind(3:4);
    TEMP1 = nanmean(HAND_TWO_V(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_V(:,ind2),2);
    
%     shuffle_V(i,1) = dtw(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
    shuffle_V(i,1) = pdist2(TEMP1',TEMP2','euclidean');%/length(TEMP1_new);
end




TEMP1_V2 = TEMP1;
TEMP2_V2 = TEMP2;

TEMP1_within1_V2 = TEMP1;
TEMP2_within1_V2 = TEMP2;



TEMP1_V = nanmean([TEMP1_V1, TEMP1_V2],2);
TEMP2_V = nanmean([TEMP2_V1, TEMP2_V2],2);

TEMP1_within1_V = nanmean([TEMP1_V1, TEMP1_V2],2);
TEMP2_within1_V = nanmean([TEMP2_V1, TEMP2_V2],2);



%%


shuffle = nanmean([shuffle_H; shuffle_V],1);
shuffle(find(shuffle>nanmean(shuffle)+2*nanstd(shuffle)))=NaN;
%  shuffle = smooth(smooth(smooth(shuffle)));




clear REAL_all MEV_all
MAV_all = [shuffle_H; shuffle_V];
REAL_all = [real_H; real_V];

VAL_p = ttest_NN(MAV_all,REAL_all); % ttest p value
VAL_R = ROC_n(MAV_all,REAL_all); % ROC value
    
% % %     TRUTH_MATRIX(inCOUNT,1) = nanmean(MAV_all);
% % %     TRUTH_MATRIX(inCOUNT,2) = nanmean(REAL_all);
% % %     TRUTH_MATRIX(inCOUNT,3) = nanmean(REAL_all)-nanmean(MAV_all);
% % % %     [~,TRUTH_MATRIX(inCOUNT,4)] = ttest2(REAL_all,MAV_all);
% % %     
% % %     TRUTH_MATRIX(inCOUNT,5) = nanmedian(MAV_all);
% % %     TRUTH_MATRIX(inCOUNT,6) = nanmedian(REAL_all);
% % %     TRUTH_MATRIX(inCOUNT,7) = nanmedian(REAL_all)-nanmedian(MAV_all);
% % % %     [TRUTH_MATRIX(inCOUNT,8),~] = ranksum(REAL_all,MAV_all);


% % % end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% % % 
% % % REAL_color = [33 178 128]/255;
% % % MAV_color = [102 102 102]/255;
% % % BIN_W = 5;
% % % F = figure;
% % % subplot(2,2,1,'box','off','tickdir','out','Linewidth',0.3,'FontSize',7)
% % % hold on;
% % % h = histogram(TRUTH_MATRIX(:,2),'BinWidth',BIN_W);
% % % h.FaceColor = REAL_color;
% % % h = histogram(TRUTH_MATRIX(:,1),'BinWidth',BIN_W);
% % % h.FaceColor = MAV_color;
% % % box off;
% % % xlabel('MRD(a.u.)');
% % % ylabel('frequency');
% % % YLIM = ylim;
% % % plot([nanmean(TRUTH_MATRIX(:,2)) nanmean(TRUTH_MATRIX(:,1))],[YLIM(2) YLIM(2)],'-k')
% % % text(nanmean(TRUTH_MATRIX(:,1))+0.25*(nanmean(TRUTH_MATRIX(:,2))-nanmean(TRUTH_MATRIX(:,1))), YLIM(2), star_n(ttest_NN(TRUTH_MATRIX(:,1),TRUTH_MATRIX(:,2))));
% % % % ylim([YLIM(1) YLIM(2)+5]);
% % % % % % % % % % % % % xlim([0 1200]);
% % % 
% % % subplot(2,2,2)
% % % hold on;
% % % text(0.5,0.5,strcat('p = ',num2str(ttest_NN(TRUTH_MATRIX(:,1),TRUTH_MATRIX(:,2)))));
% % % text(0.5,0.4,strcat('AUC = ',num2str(ROC_n(TRUTH_MATRIX(:,1),TRUTH_MATRIX(:,2)))));
% % % axis off;
% % % box off;
% % % 
% % % cd(Results_dir)
% % % filename = strcat(NOME_COMBINED,'_Motor_Learning_MAV_STATS_CENTRALTENDENCY');
% % % print(F, '-dpdf', filename, '-r400')
% % % 
% % % 
% % %     

% % % 
% % % Mot_Learn.REAL_hist = nanmean(TRUTH_MATRIX(:,2));
% % % Mot_Learn.MAV_hist = nanmean(TRUTH_MATRIX(:,1));
% % % Mot_Learn.REAL_MED = nanmean(TRUTH_MATRIX(:,6));
% % % Mot_Learn.MAV_MED = nanmean(TRUTH_MATRIX(:,5));
% % % 
% % % 
% % % Mot_Learn.VAL_p = ttest_NN(TRUTH_MATRIX(:,1),TRUTH_MATRIX(:,2));
% % % Mot_Learn.VAL_R = ROC_n(TRUTH_MATRIX(:,1),TRUTH_MATRIX(:,2));

%

REAL_color = [33 178 128]/255;
MAV_color = [102 102 102]/255;


BIN_W = 40;
F = figure;
subplot(2,2,1,'box','off','tickdir','out','Linewidth',0.3,'FontSize',7)
hold on;
h = histogram(REAL_all,'BinWidth',BIN_W);
h.FaceColor = REAL_color;
h = histogram(MAV_all,'BinWidth',BIN_W);
h.FaceColor = MAV_color;
box off;
xlabel('rms distance(a.u.)');
ylabel('frequency');
YLIM = ylim;
plot([nanmean(MAV_all) nanmean(REAL_all)],[YLIM(2) YLIM(2)],'-k')
text(nanmean(MAV_all)+0.25*(nanmean(REAL_all)-nanmean(MAV_all)), YLIM(2), star_n(VAL_p));
% ylim([YLIM(1) YLIM(2)+5]);
% % % % % % % % % % xlim([0 1200]);

subplot(2,2,2)
hold on;
text(0.5,0.5,strcat('p = ',num2str(VAL_p)));
text(0.5,0.4,strcat('AUC = ',num2str(VAL_R)));
axis off;
box off;



 filename = strcat(NOME_COMBINED,'_Association_Learning_MAV_STATS');
print(F, '-dpdf', filename, '-r400')







Mot_Learn.MAV_hist = nanmean(MAV_all);
Mot_Learn.REAL_hist = nanmean(REAL_all);
Mot_Learn.MAV_MED = nanmedian(MAV_all);
Mot_Learn.REAL_MED = nanmedian(REAL_all);
Mot_Learn.VAL_p = VAL_p;
Mot_Learn.VAL_R = VAL_R;





















TEMP1_within2_H = TEMP1;
TEMP2_within2_H = TEMP2;



F = figure();

MOT_COLOR = [209 173 21]/255;
ASS_COLOR = [33 178 127]/255;

COLOR = ASS_COLOR;

ax1 = subplot(2,3,1);
hold on;
plot(TEMP1_within1_H,'-','color',COLOR)
plot(TEMP2_within1_H,'--','color',COLOR)
axis off;
ax2 = subplot(2,3,4);
hold on;
plot(TEMP1_within1_V,'-','color',COLOR)
plot(TEMP2_within1_V,'--','color',COLOR)
axis off;
linkaxes([ax1 ax2]);
% % % 
% % % 
% % % ax1 = subplot(2,3,2);
% % % hold on;
% % % plot(TEMP1_within2_H,'-','color',COLOR)
% % % plot(TEMP2_within2_H,'--','color',COLOR)
% % % axis off;
% % % ax2 = subplot(2,3,5);
% % % hold on;
% % % plot(TEMP1_within2_V,'-','color',COLOR)
% % % plot(TEMP2_within2_V,'--','color',COLOR)
% % % axis off;
% % % linkaxes([ax1 ax2]);


COLOR = MOT_COLOR;

ax3 = subplot(2,3,3);
hold on;
plot(TEMP1_across_H,'-','color',COLOR)
plot(TEMP2_across_H,'--','color',COLOR)
axis off;
ax4 = subplot(2,3,6);
hold on;
plot(TEMP1_across_V,'-','color',COLOR)
plot(TEMP2_across_V,'--','color',COLOR)
axis off;
linkaxes([ax1 ax2 ax3 ax4]);



cd(Results_dir)
filename = strcat(NOME_COMBINED,'_skeleton');
print(F, '-dpdf', filename, '-r400');
clear F;



































REAL = nanmean(REAL_all);
SHUFFLE = nanmean(MAV_all);

Mot_Learn.REAL = REAL;
Mot_Learn.SHUFFLE = SHUFFLE;

% 
% Mot_Learn.MAV_hist = nanmean(MAV_all);
% Mot_Learn.REAL_hist = nanmean(REAL_all);
% Mot_Learn.MAV_MED = nanmedian(MAV_all);
% Mot_Learn.REAL_MED = nanmedian(REAL_all);
% Mot_Learn.VAL_p = VAL_p;
% Mot_Learn.VAL_R = VAL_R;

if strfind(NOME_COMBINED,'D')
    TASK_TYPE = 'M';   % Motor learning
else
    TASK_TYPE = 'A';   % Association learning
end 




Whole.H_mean1 = ONE.H_mean;
Whole.H_mean2 = TWO.H_mean;
Whole.V_mean1 = ONE.V_mean;
Whole.V_mean2 = TWO.V_mean;



save(POP_file,'Mot_Learn','TASK_TYPE','Whole','-append');
save(MERGE_file,'Mot_Learn','TASK_TYPE','Whole','-append');
save(ALLCELLS_file,'Mot_Learn','TASK_TYPE','Whole','-append');






%% avg waveforms -------------------




HAND_ONE = [ONE.H; ONE.V];
HAND_TWO = [TWO.H; TWO.V];

HAND_ONE_H = ONE.H;
HAND_ONE_V = ONE.V;
HAND_TWO_H = TWO.H;
HAND_TWO_V = TWO.V;









disp(' ALL DONE ');








% end