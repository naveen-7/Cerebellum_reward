%% Part of a batch process Required after CER_analysis_02
%% Analyses simple spikes in isolation
%% PRIMARILY for analysing correct and wrong trials and all their combinations for SS

% Created by NAVEEN ON 01/13/16 at CUMC
% Linearlized on 03/08/17

% function CER_analysis_03_LINEAR

% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------
% codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
% data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data';



disp('!!! CER_analysis_03_LINEAR has started running !!!');
if strcmp(ALLCELLS_file(1),'C') ALLCELLS_file(1)='E'; end
load(ALLCELLS_file);

if strcmp(POP_file(1),'C') POP_file(1)='E'; end
if strcmp(MERGE_file(1),'C') MERGE_file(1)='E'; end
if strcmp(Results_dir(1),'C') Results_dir(1)='E'; end
if strcmp(ALLCELLS_file(1),'C') ALLCELLS_file(1)='E'; end

save(POP_file,'POP_file','MERGE_file','Results_dir','Infos','-append');
save(MERGE_file,'POP_file','MERGE_file','Results_dir','Infos','-append');
save(ALLCELLS_file,'POP_file','MERGE_file','Results_dir','Infos','-append');




%% HAND CORRECTION

outcometable(:,1) = Infos(:,5);  % symbol
outcometable(:,2) = Infos(:,9); % movt
outcometable(:,3) = Infos(:,10); % reward

WHICHHAND(find(outcometable(:,2)==0 & outcometable(:,3)==1),1)=1;
WHICHHAND(find(outcometable(:,2)==1 & outcometable(:,3)==1),1)=2;
WHICHHAND(find(outcometable(:,2)==1 & outcometable(:,3)==2),1)=1;
WHICHHAND(find(outcometable(:,2)==0 & outcometable(:,3)==2),1)=2;

Infos(:,9)=WHICHHAND;







%%
%% ONE: SIMPLE SPIKES ----------------------------------------------------


disp('*************************************************');
disp('*************************************************');
disp('!!! Analysing simple spikes !!!');


%%
%%  Infos(i,10)  %Correct = 1; Wrong = 2
%%


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


SpikesS_Corr = Spikes.S;
SpikesS_Wrng = Spikes.S;
SpikesC_Corr = Spikes.C;
SpikesC_Wrng = Spikes.C;

Tgt_times_Corr = Infos(:,4);
Tgt_times_Wrng = Infos(:,4);

Mvt_times_Corr = Infos(:,11);
Mvt_times_Wrng = Infos(:,11);


for i=1:size(Infos,1)
    if Infos(i,10)==1
        SpikesS_Wrng{i,1}=[];
        SpikesC_Wrng{i,1}=[];
        Tgt_times_Wrng(i,1) = NaN;
        Mvt_times_Wrng(i,1) = NaN;
    end
    
    if Infos(i,10)==2
        SpikesS_Corr{i,1}=[];
        SpikesC_Corr{i,1}=[];
        Tgt_times_Corr(i,1) = NaN;
        Mvt_times_Corr(i,1) = NaN;
    end  
end




%% 
cd(codes_dir)

TEMP_infos = Infos; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infos(2:end,10) = Infos(1:end-1,10); % 1+N
% Infos(1:end-1,10) = Infos(2:end,10); % 2-N






% % % ANALYSE_CS
% % % 
% % % %  CER_analysis_03f_LINEAR %%%% 
% % % %  CER_analysis_03a_LINEAR %%%% 
% % % 

%    CER_analysis_03c_LINEAR  %%%% 
% %    CER_analysis_03c_LINEAR_rsm %%%% 
%    CER_analysis_03d_LINEAR  %%%% 
%    CW_DELTA_POP_crossvalidation
% %  CER_analysis_03b_LINEAR %%%% 
% % % % % 
%   CER_analysis_03g_LINEAR  %%%% 
% % % % 
% % % % %  CER_analysis_03h_LINEAR  %%%%% 
% % % % %  if ~strcmp(NOME(1),'P')
% % % % %       CER_analysis_03j_LINEAR  %%%%% 
% % % % %  end
% % % % 
% % % % %   CER_analysis_03e_LINEAR %%%%% 
%   CER_analysis_03i_LINEAR %%%% 
% % % %   CER_analysis_03k_LINEAR %%%% 
% % % 
% % % %   CER_analysis_03l_LINEAR
% % % %   CER_analysis_03m_LINEAR %%%% 
% % % CER_analysis_03m_LINEAR2 %%%% 
% % % 
% % % % CER_analysis_03n_LINEAR %%%% 
% % % 
% % % %  CER_analysis_03o_LINEAR %%% 

% CER_analysis_03pp_LINEAR

cd(codes_dir);
SS_OT

Infos = TEMP_infos; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('!!! END OF CODE 03_LINEAR !!!');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % end
