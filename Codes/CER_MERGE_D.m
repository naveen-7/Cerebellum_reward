%% To merge the rex and spike 2 files for manipulandam change task with VM change

% Written by naveen at cumc on 1/20/18

% CER_MERGE_D

close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

disp('!!!!!  CER_MERGE has started running  !!!!!')
cd('e:\NAVEEN_Work\Cerebellum\Data')
mkdir('MERGED_CELLS');
MERGE_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');
ALLCELLS_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');
POP_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Data','POP_CELLS');




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




disp('*************************************************')
disp('*************************************************')
disp('ALL FILES  WERE LOADED SUCCESSFULLY')




len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len))
MERGE_file = strcat('e:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\Data_',FileName1(1:len))
NOME = FileName1(1:len);

% getting the task type here -----------------
clear C;
[C,matches] = strsplit(NOME,{'_'});

if cell2mat(strfind(C(2),'D'))
    TASK_TYPE = 'D';
    disp('!!! TASK TYPE has been identified as D !!!');
else
    disp('!!! WARNING: Unidentified task type !!!');
    disp('!!! Please enter a task type to continue !!!');
    TASK_TYPE = upper(input('Here: ','s'));
end

clear C matches;

% ALL FILES HAVE BEEN ENTERED
% NOW CREATING A PROPER AND A SIMPLE STRUCTURE FOR FURTHER ANALYSES


cd(codes_dir);

[SPIKES.S,SHAPES.S,INFOS] = Prelim_ADAPTEDX_n(S_SPK,Trials,Trigg,'Simple');


%%%%%% FIND CHANGE

str = FileName1(1:len);
Monkey_name=upper(str(1:2));

Infos = INFOS;

cd(codes_dir);
if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    CER_analysis_02Br_CHANGE
end




disp(strcat(['----------- # trials = ',num2str(size(INFOS,1)),' ','-----------']))
disp(strcat(['----------- CHANGE   = ',num2str(CHANGE),' ','-----------']))


trials_B = (input('Enter the trial start:end for BARS : '));
trials_D = (input('Enter the trial start for DOWELS : '));
trials_L = CHANGE:size(Infos,1);

inde = randperm(length(trials_D));


clear temp;
temp = SPIKES.S;
Spikes.B = temp(trials_B);
Spikes.D = temp(trials_D);
Spikes.L = temp(trials_L);
Spikes.D = Spikes.D(inde);
Spikes.S = [Spikes.B;Spikes.D;Spikes.L];

change = length(temp(trials_B));

clear temp;
temp = SHAPES.S;
Shapes.B = temp(trials_B);
Shapes.D = temp(trials_D);
Shapes.L = temp(trials_L);
Shapes.D = Shapes.D(inde);
Shapes.S = [Shapes.B;Shapes.D;Shapes.L];

clear Infos;
INFOS1.B = INFOS(trials_B,:);
INFOS1.D = INFOS(trials_D,:);
INFOS1.L = INFOS(trials_L,:);
INFOS1.D = INFOS1.D(inde,:);
Infos = [INFOS1.B;INFOS1.D;INFOS1.L];

CHANGE = length(Spikes.B)+length(Spikes.D)+1;

disp(strcat(['----------- changed manimulandam at trial: ',num2str(change),' ','-----------']))
disp(strcat(['----------- CHANGED association at trial : ',num2str(CHANGE),' ','-----------']))



% % % %%%%% TAKING ONLY CORRECT TRIALS %%%%%
% % % clear IND
% % % IND = find(INFOS1.B(:,10)==2);
% % % Spikes.B(IND,:)=[];
% % % Shapes.B(IND,:)=[];
% % % INFOS1.B(IND,:)=[];
% % % 
% % % clear IND
% % % IND = find(INFOS1.D(:,10)==2);
% % % Spikes.D(IND,:)=[];
% % % Shapes.D(IND,:)=[];
% % % INFOS1.D(IND,:)=[];
% % % 
% % % 
% % % clear IND
% % % IND = find(Infos(:,10)==2);
% % % Spikes.S(IND,:)=[];
% % % Shapes.S(IND,:)=[];
% % % Infos(IND,:)=[];
% % % 
% % % CHANGE = length(trials_B)+1;






disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat('e:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\','RESULTS_CS_',NOME);
cd(strcat('e:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len)))
mkdir(strcat('RESULTS_CS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')



save(MERGE_file,'Infos','INFOS1','Spikes','Shapes','SPIKES','SHAPES','TASK_TYPE','NOME','CHANGE','change')
save(ALLCELLS_file,'Infos','INFOS1','Spikes','Shapes','SPIKES','SHAPES','TASK_TYPE','NOME','CHANGE','change')





%%

Monkey_name= upper(NOME(1:2));
if strcmp(Monkey_name,'BR') | strcmp(Monkey_name,'SL')
    save(MERGE_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')  
elseif strcmp(Monkey_name,'PR')
    save(MERGE_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')
    save(ALLCELLS_file,'Monkey_name','MERGE_file','ALLCELLS_file','Results_dir','POP_dir','POP_file','-append')
end










disp('ALL FILES WERE MERGED SUCCESSFULLY')
disp('!!! CELL READY FOR FURTHER ANALYSES !!!')



% end