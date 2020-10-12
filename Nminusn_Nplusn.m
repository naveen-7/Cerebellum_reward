%% Calculated delta epochs for 3-N, 2-N, 1-N, 1+N, 2+N cases

%% N-n_N+n
% % % cd(PathName1);

%% CASE

F = figure();


TEMP_infos = Infos; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Infos = Infos; CASE = '1-N'; COL = 3; % 1-N

% Infos(2:end,10) = Infos(1:end-1,10); CASE = '1+N'; COL = 4; % 1+N
% Infos(1:end-1,10) = Infos(2:end,10); CASE = '2-N'; COL = 2; % 2-N

Infos(3:end,10) = Infos(1:end-2,10); CASE = '2+N'; COL = 5; % 2+N
% Infos(1:end-2,10) = Infos(3:end,10); CASE = '3-N'; COL = 1; % 3-N



% total signal = @M(N) (200:700) + @S(N+1) (-400:200) + @M(N+1) (-300: 200)

Sigma = 50;
trials = CHANGE:CHANGE+20;

%%%%%% part1:
clear trials_corr PSTH_corr trials_wrng PSTH_wrng ind TIME
Start_time = 0; End_time = 900; TIME = Start_time:End_time;
START = 200; END = 700;
ind = find(TIME==START) : find(TIME==END); TIME = TIME(ind);
trials_corr = trials(find((Infos(trials,10))==1));
PSTH_corr = PSTH_RETURN_n(Spikes.S(trials_corr),Infos(trials_corr,11),Start_time,End_time,Sigma);
trials_wrng = trials(find((Infos(trials,10))==2));
PSTH_wrng = PSTH_RETURN_n(Spikes.S(trials_wrng),Infos(trials_wrng,11),Start_time,End_time,Sigma);

trials_OT = 1:CHANGE-1;
PSTH_OT = PSTH_RETURN_n(Spikes.S(trials_OT),Infos(trials_OT,11),Start_time,End_time,Sigma);

PSTH_corr = PSTH_corr(:,ind);
PSTH_wrng = PSTH_wrng(:,ind);
PSTH_OT = PSTH_OT(:,ind);

for i=1:length(PSTH_wrng)
    p_PART1(i) = ttest_NN(PSTH_corr(:,i),PSTH_wrng(:,i));
end


s1 = subplot(3,3,1);
hold on;
errorline_n(TIME,mean(PSTH_corr),std(PSTH_corr)/sqrt(size(PSTH_corr,1)),1,[0 0 1])
errorline_n(TIME,mean(PSTH_wrng),std(PSTH_wrng)/sqrt(size(PSTH_wrng,1)),1,[1 0 0])
errorline_n(TIME,mean(PSTH_OT),std(PSTH_OT)/sqrt(size(PSTH_OT,1)),1,[0.3 0.3 0.3])
xlim([TIME(1) TIME(end)]);
% ylim([0 100]);
ylabel('sp/s');
subplot(3,3,4)
hold on;
plot(TIME,log10(p_PART1));
xlim([TIME(1) TIME(end)]);
plot(xlim,[-1.3 -1.3],':k')
plot(xlim,[-2 -2],':k')
ylim([-10 0])
SIG = (find(p_PART1<0.05));
plot(TIME(SIG),log10(p_PART1(SIG)),'-K','linewidth',2);
xlabel('time from reward');
ylabel('log10(P)');

PSTH_C1 = PSTH_corr;
PSTH_W1 = PSTH_wrng;


%%%%%% part2:
clear trials_corr PSTH_corr trials_wrng PSTH_wrng ind TIME
Start_time = -600; End_time = 400; TIME = Start_time:End_time;
START = -400; END = 200;
ind = find(TIME==START) : find(TIME==END); TIME = TIME(ind);
trials_corr = trials(find((Infos(trials,10))==1))+1;
PSTH_corr = PSTH_RETURN_n(Spikes.S(trials_corr),Infos(trials_corr,4),Start_time,End_time,Sigma);
trials_wrng = trials(find((Infos(trials,10))==2))+1;
PSTH_wrng = PSTH_RETURN_n(Spikes.S(trials_wrng),Infos(trials_wrng,4),Start_time,End_time,Sigma);

trials_OT = 1:CHANGE-1;
PSTH_OT = PSTH_RETURN_n(Spikes.S(trials_OT),Infos(trials_OT,4),Start_time,End_time,Sigma);

PSTH_corr = PSTH_corr(:,ind);
PSTH_wrng = PSTH_wrng(:,ind);
PSTH_OT = PSTH_OT(:,ind);

for i=1:length(PSTH_wrng)
    p_PART2(i) = ttest_NN(PSTH_corr(:,i),PSTH_wrng(:,i));
end


s2 = subplot(3,3,2);
hold on;
errorline_n(TIME,mean(PSTH_corr),std(PSTH_corr)/sqrt(size(PSTH_corr,1)),1,[0 0 1])
errorline_n(TIME,mean(PSTH_wrng),std(PSTH_wrng)/sqrt(size(PSTH_wrng,1)),1,[1 0 0])
errorline_n(TIME,mean(PSTH_OT),std(PSTH_OT)/sqrt(size(PSTH_OT,1)),1,[0.3 0.3 0.3])
xlim([TIME(1) TIME(end)]);
% ylim([0 100]);
ylabel('sp/s');
subplot(3,3,5)
hold on;
plot(TIME,log10(p_PART2));
xlim([TIME(1) TIME(end)]);
plot(xlim,[-1.3 -1.3],':k')
plot(xlim,[-2 -2],':k')
ylim([-10 0])
SIG = (find(p_PART2<0.05));
plot(TIME(SIG),log10(p_PART2(SIG)),'-K','linewidth',2);
xlabel('time from symbol');
ylabel('log10(P)');

PSTH_C2 = PSTH_corr;
PSTH_W2 = PSTH_wrng;


%%%%%% part3:
clear trials_corr PSTH_corr trials_wrng PSTH_wrng
Start_time = -500; End_time = 400; TIME = Start_time:End_time;
START = -300; END = 200;
ind = find(TIME==START) : find(TIME==END); TIME = TIME(ind);
trials_corr = trials(find((Infos(trials,10))==1))+1;
PSTH_corr = PSTH_RETURN_n(Spikes.S(trials_corr),Infos(trials_corr,11),Start_time,End_time,Sigma);
trials_wrng = trials(find((Infos(trials,10))==2))+1;
PSTH_wrng = PSTH_RETURN_n(Spikes.S(trials_wrng),Infos(trials_wrng,11),Start_time,End_time,Sigma);

trials_OT = 1:CHANGE-1;
PSTH_OT = PSTH_RETURN_n(Spikes.S(trials_OT),Infos(trials_OT,11),Start_time,End_time,Sigma);


PSTH_corr = PSTH_corr(:,ind);
PSTH_wrng = PSTH_wrng(:,ind);
PSTH_OT = PSTH_OT(:,ind);

for i=1:length(PSTH_wrng)
    p_PART3(i) = ttest_NN(PSTH_corr(:,i),PSTH_wrng(:,i));
end


s3 = subplot(3,3,3);
hold on;
errorline_n(TIME,mean(PSTH_corr),std(PSTH_corr)/sqrt(size(PSTH_corr,1)),1,[0 0 1])
errorline_n(TIME,mean(PSTH_wrng),std(PSTH_wrng)/sqrt(size(PSTH_wrng,1)),1,[1 0 0])
errorline_n(TIME,mean(PSTH_OT),std(PSTH_OT)/sqrt(size(PSTH_OT,1)),1,[0.3 0.3 0.3])
xlim([TIME(1) TIME(end)]);
% ylim([0 100]);
ylabel('sp/s');
subplot(3,3,6)
hold on;
plot(TIME,log10(p_PART3));
xlim([TIME(1) TIME(end)]);
plot(xlim,[-1.3 -1.3],':k')
plot(xlim,[-2 -2],':k')
ylim([-10 0])
SIG = (find(p_PART3<0.05));
plot(TIME(SIG),log10(p_PART3(SIG)),'-K','linewidth',2);
xlabel('time from movt');
ylabel('log10(P)');


PSTH_C3 = PSTH_corr;
PSTH_W3 = PSTH_wrng;

linkaxes([s1 s2 s3],'y');


if strcmp(Results_dir(1),'C')
    Results_dir(1)='E';
end

DELTA_PROOF.USE = '1';
try
    cd(Results_dir)
catch
    TEMP = strcat(Results_dir(1:56),'_REMOVE\',Results_dir(57:end));
    cd(TEMP);
    DELTA_PROOF.USE = '0';
end
filename = strcat(NOME,'_',CASE);
print(F, '-dpdf', filename, '-r400')


%%%% temp fix 




cd('E:\NAVEEN_Work\Cerebellum\Paper\ALLCASES')
mkdir(strcat('E:\NAVEEN_Work\Cerebellum\Paper\ALLCASES\',CASE))
cd(strcat('E:\NAVEEN_Work\Cerebellum\Paper\ALLCASES\',CASE));
filename = strcat(NOME,'_',CASE);
print(F, '-dpdf', filename, '-r400')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REAL_SHUFFLE
try
PSTH_C = [PSTH_C1 PSTH_C2 PSTH_C3];
PSTH_W = [PSTH_W1 PSTH_W2 PSTH_W3];


N = 1000;

for i=1:N
    num = randperm(size(PSTH_C,1));
    num = num(1:5);
    P_C = nanmean(PSTH_C(num,:));
    
    num = randperm(size(PSTH_W,1));
    num = num(1:5);
    P_W = nanmean(PSTH_W(num,:));
    
    REAL(i) =  sqrt(nansum((P_C-P_W).^2));
    
    
    P_all = [PSTH_C;PSTH_W];
    temp = randperm(size(P_all,1));
    P_all = P_all(temp,:);
    
    P_C = nanmean(P_all(1:5,:));
    P_W = nanmean(P_all(6:10,:));
    
    SHUFF(i) =  sqrt(nansum((P_C-P_W).^2));
    
end


F = figure();
hold on;
histogram(REAL,10);
histogram(SHUFF,10);



ttest_NN(REAL,SHUFF)

if ( (lillietest(REAL)==0 & lillietest(SHUFF)==0) )
    DELTA_PROOF.REALALL(COL) = nanmean(REAL);
    DELTA_PROOF.SHUFFALL(COL) = nanmean(SHUFF);
else
    DELTA_PROOF.REALALL(COL) = nanmedian(REAL);
    DELTA_PROOF.SHUFFALL(COL) = nanmedian(SHUFF);
end
DELTA_PROOF.PALL(COL) = ttest_NN(REAL,SHUFF);

Infos = TEMP_infos; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catch
end




if strcmp(POP_file(1),'C') POP_file(1)='E'; end
if strcmp(MERGE_file(1),'C') MERGE_file(1)='E'; end
if strcmp(Results_dir(1),'C') Results_dir(1)='E'; end
try
 if strcmp(ALLCELLS_file(1),'C') ALLCELLS_file(1)='E'; end
 save(ALLCELLS_file,'DELTA_PROOF','-append');
 save(MERGE_file,'DELTA_PROOF','-append');
catch
end

try
save(POP_file,'DELTA_PROOF','-append');
catch
   save(POP_file,'DELTA_PROOF'); 
end






