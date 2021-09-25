
%% function ZERO_association_learning_manipulandam_NUEVO

% Written by naveen on 12/26/17 at cumc

% Take 5 trials before change and 5 trials after change
% take 3 random trials in bef and 3 random trials in aft and take mean
% repeat 500 times
% Randomize the trials and repeat

NUM=50;


HAND_ONE = [ONE.H; ONE.V];
HAND_TWO = [TWO.H; TWO.V];

HAND_ONE_H = ONE.H;
HAND_ONE_V = ONE.V;
HAND_TWO_H = TWO.H;
HAND_TWO_V = TWO.V;

clear real_H real_V;

% ACROSS GROUPS ------------------- H
for i=1:NUM
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_H(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_H(:,ind2),2);
    
    
    [Max,pos1] = nanmax(TEMP1); [Max,pos2] = nanmax(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);
    
    real_H(i,1) = pdist2(TEMP1_new',TEMP2_new','euclidean');
%                  dtw(TEMP1_new,TEMP2_new);
    clear TEMP1 TEMP2;
end


% ACROSS GROUPS ------------------- V
for i=1:NUM
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_V(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_V(:,ind2),2);
    
    
    [Max,pos1] = nanmin(TEMP1); [Max,pos2] = nanmin(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);
    
    real_V(i,1) = pdist2(TEMP1_new',TEMP2_new','euclidean');
%                 dtw(TEMP1_new,TEMP2_new);
    clear TEMP1 TEMP2;
end


real = nanmean([real_H; real_V],1);
real(find(real>nanmean(real)+2*nanstd(real)))=NaN;





clear shuffle_H shuffle_V;


% WITHIN GROUPS ------------------- H
for i=1:NUM
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_H(:,ind1),2);
    TEMP2 = nanmean(HAND_ONE_H(:,ind2),2);
    
    
    [Max,pos1] = nanmax(TEMP1); [Max,pos2] = nanmax(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);
    
    shuffle_H(i,1) = pdist2(TEMP1_new',TEMP2_new','euclidean');
    clear TEMP1 TEMP2;
end


TEMP1_new_H = TEMP1_new;
TEMP2_new_H = TEMP2_new;



% WITHIN GROUPS ------------------- V
for i=1:NUM/2
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_ONE_V(:,ind1),2);
    TEMP2 = nanmean(HAND_ONE_V(:,ind2),2);
    
    
    [Max,pos1] = nanmin(TEMP1); [Max,pos2] = nanmin(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);

    shuffle_V(i,1) = pdist2(TEMP1_new',TEMP2_new','euclidean');
%                 dtw(TEMP1_new,TEMP2_new);
    clear TEMP1 TEMP2;
end



TEMP1_new_V = TEMP1_new;
TEMP2_new_V = TEMP2_new;


% WITHIN GROUPS ------------------- H
for i=NUM/2+1:NUM
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_TWO_H(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_H(:,ind2),2);
    
    
    [Max,pos1] = nanmax(TEMP1); [Max,pos2] = nanmax(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);
    
    shuffle_H(i) = pdist2(TEMP1_new',TEMP2_new','euclidean');
%                 dtw(TEMP1_new,TEMP2_new);
    clear TEMP1 TEMP2;
end


% WITHIN GROUPS ------------------- V
for i=NUM/2+1:NUM
    ind1 = randperm(5); ind1=ind1(1:2);
    ind2 = randperm(5); ind2=ind2(1:2);
    TEMP1 = nanmean(HAND_TWO_V(:,ind1),2);
    TEMP2 = nanmean(HAND_TWO_V(:,ind2),2);
    
    
    [Max,pos1] = nanmin(TEMP1); [Max,pos2] = nanmin(TEMP2);
    
    if pos2>pos1
        DIFF = pos2-pos1;
        TEMP1_new = TEMP1;
        TEMP2_new = TEMP2(DIFF:end);
        L = nanmin(length(TEMP1_new),length(TEMP2_new));
        TEMP1_new=TEMP1_new(1:L);
        TEMP2_new=TEMP2_new(1:L);
    end
    
    if pos1>pos2
        DIFF = pos1-pos2;
        TEMP2_new = TEMP2;
        TEMP1_new = TEMP1(DIFF:end);
        L = nanmin(length(TEMP2_new),length(TEMP1_new));
        TEMP2_new=TEMP2_new(1:L);
        TEMP1_new=TEMP1_new(1:L);
    end
    
    TEMP1_new = TEMP1_new-TEMP1_new(1);
    TEMP2_new = TEMP2_new-TEMP2_new(1);
    
    shuffle_V(i) = pdist2(TEMP1_new',TEMP2_new','euclidean');
%                 dtw(TEMP1_new,TEMP2_new);
    clear TEMP1 TEMP2;
end




shuffle = nanmean([shuffle_H; shuffle_V],1);
shuffle(find(shuffle>nanmean(shuffle)+2*nanstd(shuffle)))=NaN;
%  shuffle = smooth(smooth(smooth(shuffle)));








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear REAL_all MEV_all
MAV_all = [shuffle_H; shuffle_V];
REAL_all = [real_H; real_V];

VAL_p = ttest_NN(MAV_all,REAL_all); % ttest p value
VAL_R = ROC_n(MAV_all,REAL_all); % ROC value

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


% cd(Results_dir)
% filename = strcat(NOME,'_Association_Learning_MAV_STATS');
% print(F, '-dpdf', filename, '-r400')


% 
% Mot_Learn.MAV_hist = nanmean(MAV_all);
% Mot_Learn.REAL_hist = nanmean(REAL_all);
% Mot_Learn.MAV_MED = nanmedian(MAV_all);
% Mot_Learn.REAL_MED = nanmedian(REAL_all);
% Mot_Learn.VAL_p = VAL_p;
% Mot_Learn.VAL_R = VAL_R;
























F = figure();

MOT_COLOR = [209 173 21]/255;
ASS_COLOR = [33 178 127]/255;

COLOR = ASS_COLOR;

ax1 = subplot(2,2,1);
hold on;
plot(TEMP1_new_H,'-','color',COLOR)
plot(TEMP2_new_H,'--','color',COLOR)
axis off;

ax2 = subplot(2,2,3);
hold on;
plot(TEMP1_new_V,'-','color',COLOR)
plot(TEMP2_new_V,'--','color',COLOR)
axis off;

linkaxes([ax1 ax2]);



% % filename = strcat(NOME_COMBINED);
% % print(F, '-dpdf', filename, '-r400');
% % clear F;




REAL = nanmean(real);
SHUFFLE = nanmean(shuffle);


figure()
subplot(3,3,5)
hold on;
bar(1,log(REAL))
bar(2,log(SHUFFLE))


figure()
subplot(3,3,5)
hold on;
plot(real)
plot(shuffle)



Mot_Learn.REAL = REAL;
Mot_Learn.SHUFFLE = SHUFFLE;


Mot_Learn.MAV_hist = nanmean(MAV_all);
Mot_Learn.REAL_hist = nanmean(REAL_all);
Mot_Learn.MAV_MED = nanmedian(MAV_all);
Mot_Learn.REAL_MED = nanmedian(REAL_all);
Mot_Learn.VAL_p = VAL_p;
Mot_Learn.VAL_R = VAL_R;

if strfind(NOME_COMBINED,'D')
    TASK_TYPE = 'M';   % Motor learning
else
    TASK_TYPE = 'A';   % Association learning
end 



save(POP_file,'Mot_Learn','TASK_TYPE','-append');
save(MERGE_file,'Mot_Learn','TASK_TYPE','-append');
save(ALLCELLS_file,'Mot_Learn','TASK_TYPE','-append');

% end