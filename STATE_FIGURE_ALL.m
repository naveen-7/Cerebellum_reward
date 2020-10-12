
% STATE_FIGURE_ALL

F = figure();

subplot(3,3,1)
hold on;
IND = 1:CHANGE;
P_MAIN_T = PSTHe_n(Spikes.S(IND),Infos(IND,4),Start_T,End_T,30,[0.5 0.5 0.5],1,1);
temp_time = Start_T:End_T;
temp_ind = find(StartLIM_T<=temp_time & temp_time<=EndLIM_T);
errorline_n(StartLIM_T:EndLIM_T,P_MAIN_T(1,temp_ind),P_MAIN_T(3,temp_ind),1,[0.5 0.5 0.5],0.3,0,1);

IND = LEARNT:size(Infos,1);
P_LEARNT_T = PSTHe_n(Spikes.S(IND),Infos(IND,4),Start_T,End_T,30,[0.5 0.5 0.5],1,1);
temp_time = Start_T:End_T;
temp_ind = find(StartLIM_T<=temp_time & temp_time<=EndLIM_T);
errorline_n(StartLIM_T:EndLIM_T,P_LEARNT_T(1,temp_ind),P_LEARNT_T(3,temp_ind),1,[0 0 0],0.3,0,1);

xlim([StartLIM_T EndLIM_T]);

MAIN_NUM = CHANGE-1;
AFTER_NUM = size(Infos,1)-LEARNT;

%%%% stats
for st=1:length(temp_ind)
    p(1,st)=ttest_n(P_MAIN_T(1,st),P_LEARNT_T(1,st),P_MAIN_T(2,st),P_LEARNT_T(2,st),MAIN_NUM,AFTER_NUM);
end
IND_p{1,1} = find(p(1,:)<0.05);
YLIM=ylim;
plot_time=Start_T:End_T;
plot(plot_time(IND_p{1,1}),[YLIM(2)*ones(length(IND_p{1,1}),1)],'.','color',[212 175 55]/255);
xlim([StartLIM_T EndLIM_T]);


subplot(3,3,2)
hold on;
IND = 1:CHANGE;
P_MAIN_M = PSTHe_n(Spikes.S(IND),Infos(IND,11),Start_M,End_M,30,[0.5 0.5 0.5],1,1);
temp_time = Start_M:End_M;
temp_ind = find(StartLIM_M<=temp_time & temp_time<=EndLIM_M);
errorline_n(StartLIM_M:EndLIM_M,P_MAIN_M(1,temp_ind),P_MAIN_M(3,temp_ind),1,[0.5 0.5 0.5],0.3,0,1);

IND = LEARNT:size(Infos,1);
P_LEARNT_M = PSTHe_n(Spikes.S(IND),Infos(IND,11),Start_M,End_M,30,[0.5 0.5 0.5],1,1);
temp_time = Start_M:End_M;
temp_ind = find(StartLIM_M<=temp_time & temp_time<=EndLIM_M);
errorline_n(StartLIM_M:EndLIM_M,P_LEARNT_M(1,temp_ind),P_LEARNT_M(3,temp_ind),1,[0 0 0],0.3,0,1);

xlim([StartLIM_M EndLIM_M]);


MAIN_NUM = CHANGE-1;
AFTER_NUM = size(Infos,1)-LEARNT;

%%%% stats
for st=1:length(temp_ind)
    p(2,st)=ttest_n(P_MAIN_M(1,st),P_LEARNT_M(1,st),P_MAIN_M(2,st),P_LEARNT_M(2,st),MAIN_NUM,AFTER_NUM);
end
IND_p{2,1} = find(p(2,:)<0.05);
YLIM=ylim;
plot_time=Start_M:End_M;
plot(plot_time(IND_p{2,1}),[YLIM(2)*ones(length(IND_p{2,1}),1)],'.','color',[212 175 55]/255);

xlim([StartLIM_M EndLIM_M]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBINs=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wind = find(Infos(CHANGE:end,10)==2)+CHANGE-1;
SPIKES =  Spikes.S(wind);
INF_T = Infos(wind,4);
INF_M = Infos(wind,11);
SIZE = size(SPIKES,1);

bins = round(linspace(1,SIZE,nBINs+1));
beginColor = [255 0 0]/255;
endColor = [139	37	0]/255;
numSteps = nBINs;
cMap = makeColorMap(beginColor, endColor, numSteps);


subplot(3,3,4)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_T = PSTH_n(SPIKES(IND),INF_T(IND),Start_T,End_T,30,cMap(i,:),0.7,0);
end
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])

subplot(3,3,5)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_M = PSTH_n(SPIKES(IND),INF_M(IND),Start_M,End_M,30,cMap(i,:),0.7,0);
end
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])


% DUP
subplot(3,3,6)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_T = PSTH_n(SPIKES(IND),INF_T(IND),Start_T,End_T,30,cMap(i,:),0.7,0);
end
subplot(3,3,9)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_M = PSTH_n(SPIKES(IND),INF_M(IND),Start_M,End_M,30,cMap(i,:),0.7,0);
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cind = find(Infos(CHANGE:end,10)==1)+CHANGE-1;
clear SPIKES INF_T INF_M
SPIKES =  Spikes.S(cind);
INF_T = Infos(cind,4);
INF_M = Infos(cind,11);
SIZE = size(SPIKES,1);

bins = round(linspace(1,SIZE,nBINs+1));
beginColor = [162	205	90]/255;
endColor = [0 0 0]/255;
numSteps = nBINs;
cMap = makeColorMap(beginColor, endColor, numSteps);


subplot(3,3,7)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_T = PSTH_n(SPIKES(IND),INF_T(IND),Start_T,End_T,30,cMap(i,:),0.7,0);
end
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])

subplot(3,3,8)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_M = PSTH_n(SPIKES(IND),INF_M(IND),Start_M,End_M,30,cMap(i,:),0.7,0);
end
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])



% DUP
subplot(3,3,6)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_T = PSTH_n(SPIKES(IND),INF_T(IND),Start_T,End_T,30,cMap(i,:),0.7,0);
end

subplot(3,3,9)
hold on;
for i=1:nBINs
    IND = bins(i):bins(i+1);
    P_MAIN_M = PSTH_n(SPIKES(IND),INF_M(IND),Start_M,End_M,30,cMap(i,:),0.7,0);
end





%%%%%%%%%% OT
oind =1:CHANGE;
clear SPIKES INF_T INF_M
SPIKES =  Spikes.S(oind);
INF_T = Infos(oind,4);
INF_M = Infos(oind,11);
SIZE = size(SPIKES,1);

subplot(3,3,6)
hold on;
P_MAIN_T = PSTH_n(SPIKES,INF_T,Start_T,End_T,30,[0 0 0],1.5,0);
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])


subplot(3,3,9)
hold on;
P_MAIN_M = PSTH_n(SPIKES,INF_M,Start_M,End_M,30,[0 0 0],1.5,0);
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])



suptitle(strcat('STATE-',NOME(1:8),'-',NOME(10:12)));


cd(Results_dir)
filename = strcat(NOME,'_STATE_CHANGE_ALLLL');
print(F, '-dpdf', filename, '-r400')




