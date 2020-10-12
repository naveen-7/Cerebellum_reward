%% CER_MERGE_K
%% To merge the rex and lick, eye files
%% written/modified by naveen at JLG aon 12/10/18




close all;
clear;
clc;

% Setup directories--------------------------------------------------------
codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('C:','NAVEEN_Work','Cerebellum','Data','RECORDED_CELLS');

disp('!!!!!  CER_MERGE_K has started running  !!!!!')
cd('C:\NAVEEN_Work\Cerebellum\Data')
mkdir('MERGED_CELLS');
MERGE_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');
ALLCELLS_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS');
POP_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Data','POP_CELLS');




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



% LOADING Spike2 - LICK FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-LICK FILE')
[FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

clear A names;
A = load(cat(2,PathName3,FileName3));
names = fieldnames(A);
LICK = A.(names{1,1});
cd(PathName1);




% LOADING Spike2 - EYEX FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-EYEX FILE')
[FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

clear A names;
A = load(cat(2,PathName3,FileName3));
names = fieldnames(A);
EYEX = A.(names{1,1});
cd(PathName1);


% LOADING Spike2 - EYEY FILE------------------------------------------------------- %% ch2
disp('*************************************************')
disp('*************************************************')
disp('LOAD Spk2-EYEY FILE')
[FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));

clear A names;
A = load(cat(2,PathName3,FileName3));
names = fieldnames(A);
EYEY = A.(names{1,1});
cd(PathName1);



% % % % LOADING Spike2 - EMPTY FILE------------------------------------------------------- %% ch2
% % % disp('*************************************************')
% % % disp('*************************************************')
% % % disp('LOAD Spk2-EMPTY FILE')
% % % [FileName3,PathName3] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
% % % disp(strcat('!!!!!','File you entered is :',FileName3,' !!!!!'));
% % %
% % % clear A names;
% % % A = load(cat(2,PathName3,FileName3));
% % % names = fieldnames(A);
% % % EMPTY = A.(names{1,1});
% % % cd(PathName1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL FILES LOADED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len))
MERGE_file = strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\Data_',FileName1(1:len))
NOME = FileName1(1:len);

% getting the task type here -----------------
clear C;
[C,matches] = strsplit(NOME,{'_'});

if cell2mat(strfind(C(2),'K'))
    TASK_TYPE = 'N';
    disp('!!! TASK TYPE has been identified as K !!!');
else
    disp('!!! WARNING: Unidentified task type !!!');
    disp('!!! Please enter a task type to continue !!!');
    TASK_TYPE = upper(input('Here: ','s'));
end

clear C matches;

% ALL FILES HAVE BEEN ENTERED
% NOW CREATING A PROPER AND A SIMPLE STRUCTURE FOR FURTHER ANALYSES



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK WAS FIGURED OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cd(codes_dir);

[LICK_TRIALS,INFOS] = Prelim_ADAPTEDK_n(LICK,Trials,Trigg);
[EYEX_TRIALS,INFOS] = Prelim_ADAPTEDK_n(EYEX,Trials,Trigg);
[EYEY_TRIALS,INFOS] = Prelim_ADAPTEDK_n(EYEY,Trials,Trigg);

% LICK_TRIALS = CER_analysis_01_LICKEYE(LICK,Trials,Trigg);
% EYEX_TRIALS = CER_analysis_01_LICKEYE(EYEX,Trials,Trigg);
% EYEY_TRIALS = CER_analysis_01_LICKEYE(EYEY,Trials,Trigg);



Infos = INFOS;

clear LICK EYE
LICK = LICK_TRIALS;
EYE.EYE_X = EYEX_TRIALS;
EYE.EYE_Y = EYEY_TRIALS;


% % % for i=1:size(EYE.EYE_X,1)
% % %     EYE.EYE_Xn{i,1} = smooth(EYE.EYE_X{i,1},0.045,'rloess');
% % %     EYE.EYE_Yn{i,1} = smooth(EYE.EYE_Y{i,1},0.045,'rloess');
% % %     i
% % % end



for i=1:size(EYE.EYE_X,1)
    
    EYE_Xvel = abs(diff(EYE.EYE_X{i,1}));
    EYE_Yvel = abs(diff(EYE.EYE_Y{i,1}));
    
    [p,pp] = findpeaks(EYE_Xvel,'Threshold',1);
    EYE_Xb = EYE.EYE_X{i,1}; EYE_Yb = EYE.EYE_Y{i,1};
    
    if length(pp)>=2
        
        if mod(length(pp),2)==0
            clear blink_pairs
            blink_pairs(1:length(pp)/2,1) = pp(1:2:length(pp));
            blink_pairs(1:length(pp)/2,2) = pp(2:2:length(pp));
            blink_pairs(1:length(pp)/2,3) = (blink_pairs(:,2)-blink_pairs(:,1));
            
            %%% removing blinks
            for j=1:length(pp)/2
                if blink_pairs(j,2)+20<length(EYE_Xb) & blink_pairs(j,1)-20>0
                    EYE_Xb(blink_pairs(j,1)-20:blink_pairs(j,2)+20)=NaN;
                    EYE_Yb(blink_pairs(j,1)-20:blink_pairs(j,2)+20)=NaN;
                elseif blink_pairs(j,2)+20>length(EYE_Xb) & blink_pairs(j,1)-20>0
                    EYE_Xb(blink_pairs(j,1)-20:length(EYE_Xb))=NaN;
                    EYE_Yb(blink_pairs(j,1)-20:length(EYE_Xb))=NaN;
                else
                    EYE_Xb(1:blink_pairs(j,2)+20)=NaN;
                    EYE_Yb(1:blink_pairs(j,2)+20)=NaN;
                end
            end
            %%% interpolating the removed values
            EYE_Xb = interpolate_NaN_n(EYE_Xb);
            EYE_Yb = interpolate_NaN_n(EYE_Yb);
            
            %%% smoothening
            EYE.EYE_Xn{i,1} = smooth(EYE_Xb,0.02,'rloess');
            EYE.EYE_Yn{i,1} = smooth(EYE_Yb,0.02,'rloess');
        else
            EYE.EYE_Xn{i,1} = [];
            EYE.EYE_Yn{i,1} = [];
        end
    else
        EYE.EYE_Xn{i,1} = EYE.EYE_X{i,1};
        EYE.EYE_Yn{i,1} = EYE.EYE_Y{i,1};
    end
    
    i
end



disp('!!! SAVING DATA !!!');
cd(MERGE_dir);

mkdir(FileName1(1:len));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len),'\','RESULTS_',NOME);
cd(strcat('C:\NAVEEN_Work\Cerebellum\Data\MERGED_CELLS\',FileName1(1:len)))
mkdir(strcat('RESULTS_',NOME))


ALLCELLS_file = strcat(ALLCELLS_dir,'\DATA_',FileName1(1:len));



POP_file = strcat(POP_dir,'\POP_',FileName1(1:len));
save(POP_file,'TASK_TYPE','NOME')


% % % save(filenm,'Infos','Spikes','-append')
if exist('LFP_raw')
    save(MERGE_file,'Infos','LICK','EYE','RAW_SIGNAL','TASK_TYPE','NOME','RAW_SIGNAL')
    save(ALLCELLS_file,'Infos','LICK','EYE','RAW_SIGNAL','TASK_TYPE','NOME','RAW_SIGNAL')
else
    save(MERGE_file,'Infos','LICK','EYE','TASK_TYPE','NOME')
    save(ALLCELLS_file,'Infos','LICK','EYE','TASK_TYPE','NOME')
end


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










disp('ALL FILES WERE MERGED SUCCESSFULLY')
disp('!!! CELL READY FOR FURTHER ANALYSES !!!')


disp('!!! Code terminated successfully !!!');
disp('!!! All data saved and ready for further processing !!!');














% end