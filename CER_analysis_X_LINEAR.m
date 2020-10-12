%% Code for analysing reward expectation


% Created by NAVEEN ON 08/10/17 at CUMC


% function CER_analysis_X_LINEAR

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
DATAfile_M = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile_M);
cd(PathName1);

Normal.Infos = Infos;
Normal.Spikes = Spikes;
Normal.NOME = NOME;



% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE "X" DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile_X = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile_X);
cd(PathName1);

Expect.Infos = Infos;
Expect.Spikes = Spikes;
Expect.NOME = NOME;



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



F = figure();


START = -1250;
END = 650;
StartLIM = START+150;
EndLIM = END-150;
temp_time = START:END;

subplot(2,2,1)
hold on;
P1 = PSTHe_n(Normal.Spikes.S,Normal.Infos(:,11),START,END,30,[1.0000    0.7255    0.0588],1,1); % Aligned to reward
P2 = PSTHe_n(Expect.Spikes.S,Expect.Infos(:,11),START,END,30,[0.6 0.6 0.6],1,1); % Aligned to reward
temp_ind = find(StartLIM<=temp_time & temp_time<=EndLIM);
errorline_n(StartLIM:EndLIM,P1(1,temp_ind),P1(3,temp_ind),1,[1.0000    0.7255    0.0588],0.3,0,1)
errorline_n(StartLIM:EndLIM,P2(1,temp_ind),P2(3,temp_ind),1,[0.6 0.6 0.6],0.3,0,1)
xlim([StartLIM EndLIM])
hold on;

MIN = nanmin([nanmin(P1(1,find(temp_time>=-1100 & temp_time<=500))-P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmin(P2(1,find(temp_time>=-1100 & temp_time<=500))-P2(3,find(temp_time>=-1100 & temp_time<=500)))])-2;

MAX = nanmax([nanmax(P1(1,find(temp_time>=-1100 & temp_time<=500))+P1(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmax(P2(1,find(temp_time>=-1100 & temp_time<=500))+P2(3,find(temp_time>=-1100 & temp_time<=500)))])+2;

ylim([MIN MAX])

plot([StartLIM StartLIM],[MIN MIN+20],'-k','linewidth',1)
plot([StartLIM StartLIM+200],[MIN MIN],'-k','linewidth',1)
axis off;

plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)

suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:11)));


%%%%%% baseline-----

START = -450;
END = 450;
StartLIM = START+150;
EndLIM = END-150;
clear temp_time
temp_time = START:END;


subplot(2,2,3)
hold on;
P3 = PSTHe_n(Normal.Spikes.S,Normal.Infos(:,3),START,END,30,[1.0000    0.7255    0.0588],1,1); % Aligned to reward
P4 = PSTHe_n(Expect.Spikes.S,Expect.Infos(:,4),START,END,30,[0.6 0.6 0.6],1,1); % Aligned to reward
temp_ind = find(StartLIM<=temp_time & temp_time<=EndLIM);
errorline_n(StartLIM:EndLIM,P3(1,temp_ind),P3(3,temp_ind),1,[1.0000    0.7255    0.0588],0.3,0,1)
errorline_n(StartLIM:EndLIM,P4(1,temp_ind),P4(3,temp_ind),1,[0.6 0.6 0.6],0.3,0,1)
xlim([StartLIM EndLIM])
hold on;

MIN = nanmin([nanmin(P3(1,find(temp_time>=-1100 & temp_time<=500))-P3(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmin(P4(1,find(temp_time>=-1100 & temp_time<=500))-P4(3,find(temp_time>=-1100 & temp_time<=500)))])-2;

MAX = nanmax([nanmax(P3(1,find(temp_time>=-1100 & temp_time<=500))+P3(3,find(temp_time>=-1100 & temp_time<=500))) ...
    nanmax(P4(1,find(temp_time>=-1100 & temp_time<=500))+P4(3,find(temp_time>=-1100 & temp_time<=500)))])+2;

ylim([MIN MAX])

plot([StartLIM StartLIM],[MIN MIN+20],'-k','linewidth',1)
plot([StartLIM StartLIM+200],[MIN MIN],'-k','linewidth',1)
axis off;

plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)






% --------------------------------------

[val,pos] = nanmax(P1(1,:));

START = -1250;
END = 650;
clear temp_time
temp_time = START:END;
Normal.PEAK = nanmean(P1(1,pos-15:pos+15));
Expect.PEAK = nanmean(P2(1,pos-50:pos+50));

START = -450;
END = 450;
clear temp_time
temp_time = START:END;
Normal.BASE = nanmean(P3(1,find(temp_time==0)-300:find(temp_time==0)-100));
Expect.BASE = nanmean(P4(1,find(temp_time==0)-300:find(temp_time==0)-100));


subplot(2,2,2)
hold on;
b1 = bar([1],[Normal.BASE]);
b1.EdgeColor=[1.0000    0.7255    0.0588];
b1.FaceColor=[1 1 1];

b2 = bar([2],[Normal.PEAK]);
b2.FaceColor=[1.0000    0.7255    0.0588];

b3 = bar([4],[Expect.BASE]);
b3.EdgeColor=[0.6 0.6 0.6];
b3.FaceColor=[1 1 1];

b4 = bar([5],[Expect.PEAK]);
b4.FaceColor=[0.6 0.6 0.6];
xlim([3-4 3+4])


cd(Results_dir)
filename = 'Reward_expectation';
print(F, '-dpdf', filename, '-r400')


REWARD.Normal_base = Normal.BASE;
REWARD.Normal_peak = Normal.PEAK;
REWARD.Expect_base = Expect.BASE;
REWARD.Expect_peak = Expect.PEAK;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear Infos Spikes;
Infos = Normal.Infos;
Spikes = Normal.Spikes;


F = figure();

subplot(4,4,1)
hist(Infos(:,17))
temp = unique(Infos(:,17));
QUANT = temp(~isnan(temp));

IND_N = find(Infos(:,17)<=5);
IND_S = union(find(Infos(:,17)==50),find(Infos(:,17)==110));
IND_L = find(Infos(:,17)>200);

COLOR_N = [0.4 0.4 0.4];
COLOR_S = [238	173	14]/255;
COLOR_L = [139	101	8]/255;




subplot(4,4,2)

START = -750;
END = 1250;
StartLIM = START+150;
EndLIM = END-150;
clear temp_time
temp_time = START:END;

clear P1 P2 P3
hold on;
P1 = PSTHe_n(Normal.Spikes.S(IND_N,1),Normal.Infos(IND_N,11),START,END,40,[0.6 0.6 0.6],1,1); 
P2 = PSTHe_n(Normal.Spikes.S(IND_S,1),Normal.Infos(IND_S,11),START,END,40,[0.6 0.6 0.6],1,1); 
P3 = PSTHe_n(Normal.Spikes.S(IND_L,1),Normal.Infos(IND_L,11),START,END,40,[0.6 0.6 0.6],1,1); 
temp_ind = find(StartLIM<=temp_time & temp_time<=EndLIM);
errorline_n(StartLIM:EndLIM,P1(1,temp_ind),P1(3,temp_ind),1,COLOR_N,0.3,0,1)
errorline_n(StartLIM:EndLIM,P2(1,temp_ind),P2(3,temp_ind),1,COLOR_S,0.3,0,1)
errorline_n(StartLIM:EndLIM,P3(1,temp_ind),P3(3,temp_ind),1,COLOR_L,0.3,0,1)
xlim([StartLIM EndLIM])
hold on;

MIN = nanmin([nanmin(P1(1,find(temp_time>=StartLIM & temp_time<=EndLIM))-P1(3,find(temp_time>=StartLIM & temp_time<=EndLIM))) ...
    nanmin(P2(1,find(temp_time>=StartLIM & temp_time<=EndLIM))-P2(3,find(temp_time>=StartLIM & temp_time<=EndLIM))) ...
    nanmin(P3(1,find(temp_time>=StartLIM & temp_time<=EndLIM))-P3(3,find(temp_time>=StartLIM & temp_time<=EndLIM)))])-2;

MAX = nanmax([nanmax(P1(1,find(temp_time>=StartLIM & temp_time<=EndLIM))+P1(3,find(temp_time>=StartLIM & temp_time<=EndLIM))) ...
    nanmax(P2(1,find(temp_time>=StartLIM & temp_time<=EndLIM))+P2(3,find(temp_time>=StartLIM & temp_time<=EndLIM))) ...
    nanmax(P3(1,find(temp_time>=StartLIM & temp_time<=EndLIM))+P3(3,find(temp_time>=StartLIM & temp_time<=EndLIM)))])-2;

ylim([MIN MAX])

plot([StartLIM StartLIM],[MIN MIN+20],'-k','linewidth',1)
plot([StartLIM StartLIM+200],[MIN MIN],'-k','linewidth',1)
axis off;

plot([0 0],[MIN MIN],'OK','MarkerFacecolor','k','markersize',4)





% --------------------------------------


clear temp_time temp_ind
temp_time = START:END;
temp_ind = find(temp_time>=400 & temp_time<=700);
rew(1) = nanmean(P1(1,temp_ind));
rew(2) = nanmean(P2(1,temp_ind));
rew(3) = nanmean(P3(1,temp_ind));


subplot(4,4,3)
hold on;
b1 = bar([1],[rew(1)]);
b1.FaceColor=COLOR_N;

b2 = bar([2],[rew(2)]);
b2.FaceColor=COLOR_S;

b3 = bar([3],[rew(3)]);
b3.FaceColor=COLOR_L;

suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:11)));


cd(Results_dir)
filename = 'Reward_AMOUNT';
print(F, '-dpdf', filename, '-r400')


REWARD.NO = rew(1);
REWARD.SMALL = rew(2);
REWARD.LARGE = rew(3);



save(POP_file,'REWARD','-append');











