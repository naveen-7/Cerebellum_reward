%% To merge the rex and spike 2 files
%% FOR CS FILES

close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('E:\NAVEEN_Work\Cerebellum\Codes\CER_codes_NEW\CS');
data_dir  = fullfile('E:\NAVEEN_Work\Cerebellum\Data\CS_CELLS\RECORDED_CS');

disp('!!!!!  CER_MERGE has started running  !!!!!')
cd('E:\NAVEEN_Work\Cerebellum\Data')
mkdir('MERGED_CELLS');
MERGE_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','CS_CELLS','MERGED_CS');
ALLCELLS_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','CS_CELLS','ALL_CS');
POP_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Data','CS_CELLS','POP_CS');




cd(data_dir);

% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE REX FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
REXfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(REXfile);
% Infos = load(REXfile);
cd(PathName1);


% LOADING Spike2 - Trigger FILE------------------------------------------------------- %% ch5
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-Trigg FILE')
[FileName2,PathName2] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName2,' !!!!!'));

clear A names;
A = load(cat(2,PathName2,FileName2));
names = fieldnames(A);
Trigg = A.(names{1,1});
cd(PathName1);



% LOADING Spike2 - SSpikes FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-S-Spikes FILE')
[FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

clear A names;
A = load(cat(2,PathName3,FileName3));
names = fieldnames(A);
S_SPK = A.(names{1,1});
cd(PathName1);



% LOADING Spike2 - CSpikes FILE-------------------------------------------------------- %% ch7
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-C-Spikes FILE')
[FileName4,PathName4] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName4,' !!!!!'));

clear A names;
A = load(cat(2,PathName4,FileName4));
names = fieldnames(A);
C_SPK = A.(names{1,1});
cd(PathName1);





%% %%%%%%%%
% 
% S_SPK.times = SPK_TIMES{1,1}/1000;
% S_SPK.values = SPK_SHAPE{1,1}/1000;
% S_SPK.codes = ones(size(S_SPK.times));
% 
% C_SPK.times = SPK_TIMES{2,1}/1000;
% C_SPK.values = SPK_SHAPE{2,1}/1000;
% C_SPK.codes = ones(size(C_SPK.times));








RAWFLAG=0;
iiter = 0;
while RAWFLAG==0
    iiter=iiter+1;
    
    if iiter==1
        disp('----------Is the RAW LFP file available?---------');
        disp('----------Please enter "y" or "n"---------');
        isRAW = upper(input('Here: ','s'));
    end
    
    if isRAW == 'Y' | isRAW == 'N'
        RAWFLAG=1;
    else
        disp('----------Please enter a valid input--------');
        disp('---------- Please enter "y" or "n" ---------');
        isRAW = upper(input('Here: ','s'));
    end
end



if isRAW == 'Y'
    %     LOADING Spike2 - Raw FILE------------------------------------------------------- %% ch1
    disp('*************************************************')
    disp('*************************************************')
    disp('LOAD Spk2-LFP FILE')
    [FileName5,PathName5] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
    disp(strcat('!!!!!','File you entered is :',FileName5,' !!!!!'));
    
    clear A names;
    A = load(cat(2,PathName5,FileName5));
    names = fieldnames(A);
    CH1 = A.(names{1,1});
    cd(PathName1);
    
    
elseif isRAW =='N'
    disp('----------No RAW LFP file was inputed---------');
    disp('>>>>>>>>>>>>      Proceeding     >>>>>>>>>>>>>');
end



disp('*************************************************')
disp('*************************************************')
disp('ALL FILES  WERE LOADED SUCCESSFULLY')




len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len))
MERGE_file = strcat(MERGE_dir,'\',FileName1(1:len),'\Data_',FileName1(1:len));
NOME = FileName1(1:len);

% getting the task type here -----------------
clear C;
[C,matches] = strsplit(NOME,{'_'});

if cell2mat(strfind(C(2),'N'))
    TASK_TYPE = 'N';
    disp('!!! TASK TYPE has been identified as N !!!');
elseif cell2mat(strfind(C(2),'R'))
    TASK_TYPE = 'R';
    disp('!!! TASK TYPE has been identified as R !!!');
elseif cell2mat(strfind(C(2),'W'))
    TASK_TYPE = 'W';
    disp('!!! TASK TYPE has been identified as W !!!');
elseif cell2mat(strfind(C(2),'X'))
    TASK_TYPE = 'X';
    disp('!!! TASK TYPE has been identified as X !!!');
else
    disp('!!! WARNING: Unidentified task type !!!');
    disp('!!! Please enter a task type to continue !!!');
    TASK_TYPE = upper(input('Here: ','s'));
end

clear C matches;

% ALL FILES HAVE BEEN ENTERED
% NOW CREATING A PROPER AND A SIMPLE STRUCTURE FOR FURTHER ANALYSES


cd(codes_dir);

[Spikes.S,Shapes.S,INFOS] = Prelim_ADAPTED_n(S_SPK,Trials,Trigg,'Simple');
[Spikes.C,Shapes.C,Infos] = Prelim_ADAPTED_n(C_SPK,Trials,Trigg,'Complex');



if isRAW == 'Y'
    CER_analysis_01_RAW;
    RAW_SIGNAL;
end


Infos = INFOS;


disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat(MERGE_dir,'\',FileName1(1:len),'\','RESULTS_CS_',NOME);
cd(strcat(MERGE_dir,'\',FileName1(1:len)))
mkdir(strcat('RESULTS_CS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')


% % % save(filenm,'Infos','Spikes','-append')
if exist('LFP_raw')
    save(MERGE_file,'Infos','Spikes','Shapes','RAW_SIGNAL','TASK_TYPE','NOME','RAW_SIGNAL')
    save(ALLCELLS_file,'Infos','Spikes','Shapes','RAW_SIGNAL','TASK_TYPE','NOME','RAW_SIGNAL')
else
    save(MERGE_file,'Infos','Spikes','Shapes','TASK_TYPE','NOME')
    save(ALLCELLS_file,'Infos','Spikes','Shapes','TASK_TYPE','NOME')
end

% 
%  save(MERGE_file,'Infos','Spikes','Shapes','TASK_TYPE','NOME','-append')
%     save(ALLCELLS_file,'Infos','Spikes','Shapes','TASK_TYPE','-append')
%     
    
    

%%

% Getting basic cell info ---------

str = FileName1(1:len);
Monkey_name=upper(str(1:2));

cd(codes_dir);
if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    CER_analysis_02Br_CHANGE
elseif strcmp(Monkey_name,'PR')
    CER_analysis_02_CHANGE
end

% close all;

CHANGE



if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    save(MERGE_file,'Monkey_name','MERGE_file','CHANGE','LEARNT','ALL_CHANGES','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','CHANGE','LEARNT','ALL_CHANGES','Results_dir','POP_dir','POP_file','-append')  
elseif strcmp(Monkey_name,'PR')
    save(MERGE_file,'Monkey_name','MERGE_file','CHANGE','LEARNT','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','CHANGE','LEARNT','Results_dir','POP_dir','POP_file','-append')
end









%% CLEAN UP CS ------------------------------------------



Proove_CS_RANDOM_n

% % % TRY TO FIND SPURIOUS CS and elimanate them
% %
% % [S_SPK_aligned First_SPK] = Align_CS_n(SPIKES_Main.S,SPIKES_Main.C,Start,End);
% %
% % nanmean(First_SPK)-(3*(nanstd(First_SPK)))
% % nanmedian(First_SPK)


disp('ALL FILES WERE MERGED SUCCESSFULLY')
disp('!!! CELL READY FOR FURTHER ANALYSES !!!')


isrun = upper(input('Do you want to run PRELIM_ALL_n?  ','s'));
if strcmp(isrun,'Y')
    PRELIM_ALL_n
end



disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp(strcat('!!! Is this a good learning?'));
LEARNING =upper(input('Enter here: ','s'));


save(MERGE_file,'LEARNING','-append')
save(ALLCELLS_file,'LEARNING','-append')
save(POP_file,'LEARNING','-append')


disp('!!! Code terminated successfully !!!');
disp('!!! All data saved and ready for further processing !!!');














% end