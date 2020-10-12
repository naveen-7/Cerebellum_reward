

function HAND_MERGE_n




clc;
clear ;
close all;


codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES';


disp('!!!!!  HAND_MERGE_n has started running  !!!!!')
cd(data_dir)

% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('load the FIRST data file')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATA_file = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATA_file)
cd(PathName1);



NOME = FileName1(1:12);


cd('C:\NAVEEN_Work\Cerebellum\Data\videos\MERGED')






%% 
%%%%%%%%%% doing first %%%%%%%%%%%%%%%%

hand = HAND_new;
% hand = HAND_chosen;
% [H1,H2] = alignsignals(hand.H{1,1},hand.H{1,2});
% [X2,H3] = alignsignals(hand.H{1,1},hand.H{1,3});
% [X3,H4] = alignsignals(hand.H{1,1},hand.H{1,4});
% [X4,H5] = alignsignals(hand.H{1,1},hand.H{1,5});

H11 = hand.H{1,1};
H21 = hand.H{1,2};
H31 = hand.H{1,3};
H41 = hand.H{1,4};
H51 = hand.H{1,5};


% H1 = H1(find(H1>0,1):end);
% H2 = H2(find(H2>0,1):end);
% H3 = H3(find(H3>0,1):end);
% H4 = H4(find(H4>0,1):end);
% H5 = H5(find(H5>0,1):end);

MIN1 = nanmin([length(H11) length(H21) length(H31) length(H41) length(H51)]);
% H11 = H1(1:MIN); H21 = H2(1:MIN); H31 = H3(1:MIN); H41 = H4(1:MIN); H51 = H5(1:MIN);
% H_MEAN = nanmean([H1 H2 H3 H4 H5],2);


% [V1,V2] = alignsignals(hand.V{1,1},hand.V{1,2});
% [X2,V3] = alignsignals(hand.V{1,1},hand.V{1,3});
% [X3,V4] = alignsignals(hand.V{1,1},hand.V{1,4});
% [X4,V5] = alignsignals(hand.V{1,1},hand.V{1,5});

V11 = hand.V{1,1};
V21 = hand.V{1,2};
V31 = hand.V{1,3};
V41 = hand.V{1,4};
V51 = hand.V{1,5};

% 
% V1 = V1(find(V1>0,1):end);
% V2 = V2(find(V2>0,1):end);
% V3 = V3(find(V3>0,1):end);
% V4 = V4(find(V4>0,1):end);
% V5 = V5(find(V5>0,1):end);

MIN2 = nanmin([length(V11) length(H21) length(V31) length(V41) length(V51)]);
% V11 = V1(1:MIN); V21 = V2(1:MIN); V31 = V3(1:MIN); V41 = V4(1:MIN); V51 = V5(1:MIN);
% V_MEAN = nanmean([V1 V2 V3 V4 V5],2);

NOME_ONE = NOME;

disp('!!!!!  File one successfully added  !!!!!')





%% 
%%%%%%%%%% doing next %%%%%%%%%%%%%%%%

clc;
clearvars -except F ONE NOME_ONE codes_dir data_dir ax1 ax2 MIN1 MIN2 ...
           H11 H21 H31 H41 H51 V11 V21 V31 V41 V51;

cd(data_dir)

% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('load the SECOND data file')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATA_file = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATA_file)
cd(PathName1);



NOME = FileName1(1:12);
cd('C:\NAVEEN_Work\Cerebellum\Data\videos\MERGED')


hand = HAND_new;
% [X3,H4] = alignsignals(hand.H{1,1},hand.H{1,4});
% hand = HAND_chosen;
% [H1,H2] = alignsignals(hand.H{1,1},hand.H{1,2});

% [X2,H3] = alignsignals(hand.H{1,1},hand.H{1,3});
% [X4,H5] = alignsignals(hand.H{1,1},hand.H{1,5});
H12 = hand.H{1,1};
H22 = hand.H{1,2};
H32 = hand.H{1,3};
H42 = hand.H{1,4};
H52 = hand.H{1,5};


% H1 = H1(find(H1>0,1):end);
% H2 = H2(find(H2>0,1):end);
% H3 = H3(find(H3>0,1):end);
% H4 = H4(find(H4>0,1):end);
% H5 = H5(find(H5>0,1):end);

MIN3 = nanmin([length(H12) length(H22) length(H32) length(H42) length(H52)]);
% H12 = H1(1:MIN); H22 = H2(1:MIN); H32 = H3(1:MIN); H42 = H4(1:MIN); H52 = H5(1:MIN);
% H_MEAN = nanmean([H1 H2 H3 H4 H5],2);


% [V1,V2] = alignsignals(hand.V{1,1},hand.V{1,2});
% [X2,V3] = alignsignals(hand.V{1,1},hand.V{1,3});
% [X3,V4] = alignsignals(hand.V{1,1},hand.V{1,4});
% [X4,V5] = alignsignals(hand.V{1,1},hand.V{1,5});

V12 = hand.V{1,1};
V22 = hand.V{1,2};
V32 = hand.V{1,3};
V42 = hand.V{1,4};
V52 = hand.V{1,5};


% V1 = V1(find(V1>0,1):end);
% V2 = V2(find(V2>0,1):end);
% V3 = V3(find(V3>0,1):end);
% V4 = V4(find(V4>0,1):end);
% V5 = V5(find(V5>0,1):end);

MIN4 = nanmin([length(V12) length(H22) length(V32) length(V42) length(V52)]);
% V12 = V1(1:MIN); V22 = V2(1:MIN); V32 = V3(1:MIN); V42 = V4(1:MIN); V52 = V5(1:MIN);
% V_MEAN = nanmean([V1 V2 V3 V4 V5],2);


NOME_TWO = NOME;

disp('!!!!!  File two successfully added  !!!!!')



%%% CUT ALL ---------------------------------

MIN = nanmin([MIN1 MIN2 MIN3 MIN4]);

H11 = H11(1:MIN);   V11 = V11(1:MIN);
H21 = H21(1:MIN);   V21 = V21(1:MIN);
H31 = H31(1:MIN);   V31 = V31(1:MIN);
H41 = H41(1:MIN);   V41 = V41(1:MIN);
H51 = H51(1:MIN);   V51 = V51(1:MIN);

H12 = H12(1:MIN);   V12 = V12(1:MIN);
H22 = H22(1:MIN);   V22 = V22(1:MIN);
H32 = H32(1:MIN);   V32 = V32(1:MIN);
H42 = H42(1:MIN);   V42 = V42(1:MIN);
H52 = H52(1:MIN);   V52 = V52(1:MIN);





ONE.H = [H11 H21 H31 H41 H51];
ONE.V = [V11 V21 V31 V41 V51];

H_MEAN1 = nanmean([H11 H21 H31 H41 H51],2);
V_MEAN1 = nanmean([V11 V21 V31 V41 V51],2);
ONE.H_mean = H_MEAN1;
ONE.V_mean = V_MEAN1;


TWO.H = [H12 H22 H32 H42 H52];
TWO.V = [V12 V22 V32 V42 V52];

H_MEAN2 = nanmean([H12 H22 H32 H42 H52],2);
V_MEAN2 = nanmean([V12 V22 V32 V42 V52],2);
TWO.H_mean = H_MEAN2;
TWO.V_mean = V_MEAN2;




%%%%% FIGURE




F = figure();
ax1 = subplot(2,2,1);
hold on;
plot(H11); plot(H21); plot(H31); plot(H41); plot(H51);
plot(H_MEAN1,'-K','linewidth',1.5)
xlim([1 length(H11)])
title(NOME_ONE)

ax2 = subplot(2,2,3);
hold on;
plot(V11); plot(V21); plot(V31); plot(V41); plot(V51);
plot(V_MEAN1,'-K','linewidth',1.5)
xlim([1 length(V11)])

linkaxes([ax1,ax2],'xy')



ax3 = subplot(2,2,2);
hold on;
plot(H12); plot(H22); plot(H32); plot(H42); plot(H52);
plot(H_MEAN2,'-K','linewidth',1.5)
xlim([1 length(H11)])
title(NOME_TWO)

ax4 = subplot(2,2,4);
hold on;
plot(V12); plot(V22); plot(V32); plot(V42); plot(V52);
plot(V_MEAN2,'-K','linewidth',1.5)
xlim([1 length(V11)])

linkaxes([ax1,ax2,ax3,ax4],'xy')


CORR_H = corr(H_MEAN1,H_MEAN2)
CORR_V = corr(V_MEAN1,V_MEAN2)












%% SAVING STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except F ONE TWO NOME_ONE NOME_TWO codes_dir data_dir;




NOME_COMBINED = strcat(NOME_ONE,'_',NOME_TWO);


mkdir(strcat(NOME_COMBINED));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\MERGED\',NOME_COMBINED);
cd(Results_dir)

filename = NOME_COMBINED;
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');

clear F;


MERGE_dir = ('C:\NAVEEN_Work\Cerebellum\Data\videos\MERGED');
ALLCELLS_dir = ('C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES');
POP_dir = ('C:\NAVEEN_Work\Cerebellum\Data\videos\POPULATION');



MERGE_file = strcat(MERGE_dir,'\Data_',NOME_COMBINED);
ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',NOME_COMBINED);
POP_file = strcat(POP_dir,'\POP_',NOME_COMBINED);

save(MERGE_file,'ONE','TWO','NOME_COMBINED','POP_file','MERGE_file','ALLCELLS_file')
save(ALLCELLS_file,'ONE','TWO','NOME_COMBINED','POP_file','MERGE_file','ALLCELLS_file')
save(POP_file,'ONE','TWO','NOME_COMBINED','POP_file','MERGE_file','ALLCELLS_file')



disp('!!!!!  ALL CONVERSIONS DONE  !!!!!')





end
