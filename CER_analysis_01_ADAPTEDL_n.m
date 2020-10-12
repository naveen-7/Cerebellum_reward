%% Takes data from REX and Spike2 and merges them and makes the data ready
%% for further analysis

% Created by Anna and Adam
% Modified by Naveen on 11/30/15 at CUMC
% Adapted by Naveen on 12/29/16 at CUMC
% Xd by Naveen on 08/11/17 at CUMC
% Ld by Naveen on 08/23/17 at CUMC

% function CER_analysis_01_ADAPTEDL_n


% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------

% cd(data_dir)
disp('!!!!!  CER_analysis_01_ADAPTEDL_n has started running  !!!!!')


% CODES ----------------------------------------------------

START_ = 800;
FIXON_ = 1003;
TARG_ = 1007;
CORR_ = 1012;
WRONG_ = 1013;
EOT_ = 1014; 
STARTSPKREC_ = 1095;
ENDSPKREC_ = 1096;
END_ = 801;




%%
%% ONE: ALIGN REX AND SPIKE2
disp('*************************************************')
disp('*************************************************')
disp('ALIGNING REX AND SPIKE2')

data={};        % final data structure
% column 1: unsorted spike times
% column 2: event times
% column 3: sorted spike times
% column 4: rex sorted
Rexdiff=[];     % time differences between 1095 and 1096 for each trial from REX
Spikediff=[];   % time differences between 1095 and 1096 for each trial from Spike2
diffdiff=[];    % the difference between the two previously mentioned differences



% DO NOT TAKE FIRST TRIAL
% find time of first 1001 eCode to calculate start time of first trial
startTimes=[NaN];
endTimes=[NaN];
loc95=NaN;
for i=size(Trials(1).Events,2):-1:1
    if Trials(1).Events(i).Code == START_
        loc95=i+1;
        startTimes=[double(Trials(1).Events(i).Time)];
        endTimes=[double(Trials(1).Events(size(Trials(1).Events,2)).Time)];
        Rexdiff=[double(Trials(1).Events(size(Trials(1).Events,2)-1).Time)-double(Trials(1).Events(i+1).Time)];
        break;
    end
end

for i=2:size(Trials,2)
    startTimes=[startTimes; double(Trials(i).Events(1).Time)];
    endTimes=[endTimes; double(Trials(i).Events(size(Trials(i).Events,2)).Time)];
    Rexdiff=[Rexdiff; double(Trials(i).Events(size(Trials(i).Events,2)-1).Time)-double(Trials(i).Events(2).Time)];
end

for i=2:2:size(Trigg.times,1)
    Spikediff=[Spikediff; Trigg.times(i)-Trigg.times(i-1)];
end
Spikediff=Spikediff*1000;



% check if we have collected the wrong intervals from Spike
count = 0;
start3=0;
for i=1:size(Spikediff,1)
    if Spikediff(i) > 779 && Spikediff(i) < 801
        count = count+1;
    else
        count = 0;
    end
    if count >= 5
        Spikediff=[];
        for i=3:2:size(Trigg.times,1)
            Spikediff=[Spikediff; Trigg.times(i)-Trigg.times(i-1)];
        end
        Spikediff=Spikediff*1000;
        start3=1;
        break;
    end
end




sBig=[];
if size(Spikediff,1) > size(Rexdiff,1)
    sBig=1;
    A=Spikediff;
    B=Rexdiff;
else
    sBig=0;
    B=Spikediff;
    A=Rexdiff;
end

%find trial offset between REX and Spike
exit = 0;
match=0;
offset=0; % how far the shorter array of A/B is shifted to match the longer array
while match==0 && exit == 0
    diffdiff=A(offset+1:offset+size(B,1))-B;
    for i=1:20
        if diffdiff(i)>1
            offset=offset+1;
            if offset+size(B,1)>size(A,1)
                exit = 1;
            end
            break;
        end
        if i==20
            match=1;
        end
    end
end
% 
% if match==0
%     offset=-1;
%     while match==0
%         diffdiff=A(1:size(B,1)+offset)-B(1-offset:size(B,1));
%         for i=1:20
%             if diffdiff(i)>1
%                 offset=offset-1;
%                 break;
%             end
%             if i==20
%                 match=1;
%             end
%         end
%     end
% end


% THESE LINES ARE FOR MANUAL CORRECTION
% offset=10;
% diffdiff=A(offset+1:offset+size(B,1))-B;

timeOffset=[]; % Rex-1000*Spike.  Time from 1000*Spike to Rex
if sBig && offset < 0
    timeOffset=double(Trials(1-offset).Events(2).Time)-Trigg.times(1+start3)*1000;
elseif sBig && offset >= 0
    timeOffset=double(Trials(1).Events(loc95).Time)-Trigg.times(2*offset+1+start3)*1000;
elseif ~sBig && offset >= 0
    timeOffset=double(Trials(1+offset).Events(2).Time)-Trigg.times(1+start3)*1000;
elseif ~sBig && offset < 0
    timeOffset=double(Trials(1).Events(loc95).Time)-Trigg.times(2*(-offset)+1+start3)*1000;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LICK SORTING
disp('*************************************************')
disp('*************************************************')



clear DATA_LICK
temp_LICK = LICK_DATA_IMPORT.data;
DATA_LICK = temp_LICK(find(temp_LICK(:,2)==1),1);

DATA_LICK_fin = 1./diff(DATA_LICK);

clear temp_ind
for i=2:length(DATA_LICK)
    temp_ind(i-1,1) = find(temp_LICK(:,1)==DATA_LICK(i));
end

TEMP = NaN(size(temp_LICK,1),1);
TEMP(temp_ind)=DATA_LICK_fin;

x = TEMP;
nanx = isnan(x);
t    = 1:numel(x);
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));

x(isnan(x))=0;

CHAN.values = x;
CHAN.times = temp_LICK(:,1);

    
%%

trial=1;
spikes=[];
all_shape=[];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
types=1;
sorted=cell(types,1);
initial_shape=cell(types,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%
LICK={}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;

while i<size(CHAN.times,1)
    if CHAN.times(i)*1000+timeOffset<startTimes(1)
        i=i+1;
    else
        break;
    end
end


while i<size(CHAN.times,1)
     trial
    if (trial+1)>size(startTimes,1) || CHAN.times(i)*1000+timeOffset<startTimes(trial+1)
%         spikes=[spikes;CHAN.times(i)*1+timeOffset-startTimes(trial)];
        all_shape=[all_shape;CHAN.values(i,:)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        LICK{trial,1}=all_shape;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trial=trial+1;
        all_shape=[];                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        continue;
    end
    i=i+1;
end

LICK{trial,1}=all_shape;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  THREE Getting all the events
disp('*************************************************')
disp('*************************************************')
disp('GETTING ALL THE EVENTS')

for k=1:size(data,1)
    datcode=data{k,2};
    sp=data{k,3};
    %     datsp=round(sp{1})
    %     [the_struct(xx).Events.Code; the_struct(xx).Events.Time]
    dat(k).Events.Code=[datcode.Code];
    dat(k).Events.Time=double([datcode.Time]);

    for y=1:length(data{k,3})
        
        dat(k).Units(y).Code=y;
        datsp=round(sp{y});
        dat(k).Units(y).Times=double(round(datsp))';
    end
    eventpack= [double([datcode.Code])' double([datcode.Time])'];
end








% filenm=[filename '_SSspk']
leng = length(Trials); % the # trials without trial 1.
for xx = 1:leng
    

% xx
    event_packet = [double([Trials(xx).Events.Code]); double([Trials(xx).Events.Time])]';
    num_events = size(event_packet,1);
    if  num_events <=40 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for yy = 2:num_events
            Trials(xx).Events(yy).Time=(Trials(xx).Events(yy).Time-Trials(xx).Events(1).Time);
        end

        event_packets=[[Trials(xx).Events.Code]' [Trials(xx).Events.Time]' ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_units = length(Trials(xx).Units);%.Times);
        %     if Trials(xx).Units.Times(1)>10000
        for yy = 1:num_units
            Trials(xx).Units(yy).Times=(Trials(xx).Units(yy).Times-Trials(xx).Events(1).Time);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Put in array number
        if ~isempty(event_packet(find((event_packet(:,1) >= 2000) & (event_packet(:,1) < 3000)),1) - 2000)
            Trials(xx).Array = event_packet(find((event_packet(:,1) >= 2000) & (event_packet(:,1) < 3000)),1) - 2000;
        else
            Trials(xx).Array=NaN;
        end

        if ~isempty (Trials(xx).Events(find(event_packet(:,1) == 1020)))
            Trials(xx).barrelease=1;
        end
        % quewto loop "if" ho dovuto aggiungerlo perche'
        % molti trials non c'era il code 1007?????????
        % 16 aprile 02
        if ~isempty (Trials(xx).Events(find(event_packet(:,1) == 1007)))
            rf_on_time = double(event_packet(find(event_packet(:,1)==1007),2));
            if size(rf_on_time,1)>1
                rf_on_time=rf_on_time(1);
            end
            %             if rf_on_time > 0

            % Fix y-eye position by shift_amount
            % This ASSUMES that the eye postion (eyev) is signal # 2.
            if ~isempty(Trials(xx).Signals)
                temp = [double(Trials(xx).Signals(2).Signal)]';
                temp = temp;% + shift_amount;
                for yy = 1:(length(temp))
                    Trials(xx).Signals(2).Signal(yy) = temp(yy);
                end
            end

            % Now find RF object information for each trial
            first_patcd = min(find((event_packet(:,1)>3000)& (event_packet(:,1)<4000) &(event_packet(:,2)==rf_on_time))); % The 3000 code gives the # of objects following
            num_objs = double(event_packet(first_patcd,1)) - 3000;
            last_patcd = max(find((event_packet(:,2) == rf_on_time)&(event_packet(:,1) > 5000))); % this is to double check

            if (mod(((last_patcd + 1) - (first_patcd + 1)),6)>0)
                %err_text = ['In Trial ' num2str(xx) 'There are ' num2str((last_patcd + 1) - first_patcd) ' objects.'];
                %errordlg('There are not an even number of codes','Warning');
            else
                if ((((last_patcd + 1) - (first_patcd + 1))/6)~= num_objs)
                    %  errordlg('There are a different number of rf_objects to what is specificed','Warning');
                end
            end

            if isempty(first_patcd)

                Trials(xx).obj.pattern = NaN;
                Trials(xx).obj.xpos = NaN;
                Trials(xx).obj.ypos = NaN;
                Trials(xx).obj.red = NaN;
                Trials(xx).obj.green = NaN;
                Trials(xx).obj.blue = NaN;
 
            end
            
            
            counter = first_patcd + 1; % because first code is 3000 (number of objects)
            for yy = 1:num_objs
                Trials(xx).obj(yy).pattern = event_packet(counter,1) - 7000;
                counter = counter + 1;
                Trials(xx).obj(yy).xpos = (event_packet(counter,1) - 6000)/10;
                counter = counter + 1;
                Trials(xx).obj(yy).ypos = (event_packet(counter,1) - 6500)/10;
                counter = counter + 1;
                Trials(xx).obj(yy).red = event_packet(counter,1) - 5000;
                counter = counter + 1;
                Trials(xx).obj(yy).green = event_packet(counter,1) - 5300;
                counter = counter + 1;
                Trials(xx).obj(yy).blue = event_packet(counter,1) - 5600;
                counter = counter + 1;
            end
            if ~isempty(Trials(xx).obj)%(1).pattern);

                Trials(xx).pattern=Trials(xx).obj(1).pattern;
            else 
                            Trials(xx).pattern=NaN;
            end
            Trials(xx).rf_on_time = double(event_packet(find(event_packet(:,1)==1007),2));
            if size(Trials(xx).rf_on_time,1)>1
                 Trials(xx).rf_on_time= Trials(xx).rf_on_time(1);
            end
            if ~isempty(event_packet(find(event_packet(:,1)==1010),2))
                Trials(xx).rf_off_time = double(event_packet(find(event_packet(:,1)==1010),2));
             else
                 Trials(xx).rf_off_time = NaN;
            end
            
            
            if ~isempty(event_packet(find(event_packet(:,1)==1093),2))
                Trials(xx).abort=double(event_packet(find(event_packet(:,1)==1093),2));
             else
                 Trials(xx).abort=NaN;
            end


            if ~isempty(event_packet(find(event_packet(:,1)==1003),2))
                Trials(xx).FP_on=double(event_packet(find(event_packet(:,1)==1003),2));
            else
                Trials(xx).FP_on=NaN;
            end

            if ~isempty(event_packet(find(event_packet(:,1)==1012),2))
                Trials(xx).correct = 1;
                Trials(xx).reward= double(event_packet(find(event_packet(:,1)==1012),2));%-double([Trials(xx).Events(1).Time]);
                Trials(xx).MRT=Trials(xx).reward-Trials(xx).rf_on_time;
            else
                Trials(xx).correct = NaN;
                Trials(xx).reward= NaN;
                Trials(xx).MRT=NaN;
            end
            
            if ~isempty(event_packet(find(event_packet(:,1)==1013),2)) ...
                & isempty(event_packet(find(event_packet(:,1)==1093),2)) ... %%%%%% CHANGED BY NAVEEN ON 2/10
                 & isempty(event_packet(find(event_packet(:,1)==1090),2))    %%% WRONG = 1013 without 1093 and 1090
                Trials(xx).correct = 0;
                Trials(xx).reward=NaN;
                error_time=double(event_packet(find(event_packet(:,1)==1013),2));
                Trials(xx).MRT=error_time-Trials(xx).rf_on_time;
            end

        else
            Trials(xx).pattern=NaN;
            Trials(xx).rf_on_time =NaN;
            Trials(xx).rf_off_time = NaN;
            Trials(xx).abort=NaN;
            Trials(xx).FP_on=NaN;
            Trials(xx).correct = NaN;
            Trials(xx).reward= NaN;
            Trials(xx).MRT=NaN;
            Trials(xx).obj.pattern = NaN;
            Trials(xx).obj.xpos = NaN;
            Trials(xx).obj.ypos = NaN;
            Trials(xx).obj.red = NaN;
            Trials(xx).obj.green = NaN;
            Trials(xx).obj.blue = NaN;
            Trials(xx).Array=NaN;

        end
    else
        Trials(xx).pattern=NaN;
        Trials(xx).rf_on_time =NaN;
        Trials(xx).rf_off_time = NaN;
        Trials(xx).abort=NaN;
        Trials(xx).FP_on=NaN;
        Trials(xx).correct = NaN;
        Trials(xx).reward= NaN;
        Trials(xx).MRT=NaN;
        Trials(xx).obj.pattern = NaN;
        Trials(xx).obj.xpos = NaN;
        Trials(xx).obj.ypos = NaN;
        Trials(xx).obj.red = NaN;
        Trials(xx).obj.green = NaN;
        Trials(xx).obj.blue = NaN;
        Trials(xx).Array=NaN;

    end

end














%% NOW MAKE ALL SORT AND SAVE
%% find NAN array
ar=[Trials.Array]';
mrt=[Trials.MRT]';
ptr=[Trials.pattern]';

thistr=(~isnan(ar(:,1)) & ~isnan(mrt(:,1)) & ~isnan(ptr(:,1)));
Trials=Trials(thistr);

LICK=LICK(thistr,:);
% dat=dat(thistr);


%% in some trials there is not the Units Code; take out those trials
wr=1;
for u=1:length(Trials);
    if ~isempty(Trials(u).Units);
        Trials(wr)=Trials(u);
        wr=wr+1;
    end
    allTrials(u,1)=Trials(u);
end



% First find only completed correct trials
good_trials = search_for(Trials,1012);
if ~isempty(good_trials)
    good_trials = good_trials(:,1);
end

rfon_trials = search_for(Trials,1007);
bad_trials=search_for(Trials,1013);
bar_release=search_for(Trials,1020);
abort_trial=search_for(Trials,1090);
fix_break=search_for(Trials,1091);

Wrong_tri=[];
Good_tri=[];
if ~isempty(bad_trials)
    thiswrong=intersect(rfon_trials(:,1),bad_trials(:,1));
    tr=1;
    for c=1:length(thiswrong)
        if ~isempty(Trials(thiswrong(c)).MRT)
            wrong_bar(tr,1)=thiswrong(c);
            tr=tr+1;
        end
    end

    %% in some trials MRT is empty even if the trialsd started
    %% clean the trials
    wrong_mrt= [Trials(wrong_bar).MRT]';
    bad_array=[Trials(wrong_bar).Array]';
    bad_trials_patterns=[Trials(wrong_bar).pattern]';
    Wrong_tri=[wrong_bar wrong_mrt bad_array bad_trials_patterns]; %%%%%%
    Wrong_tri(:,end+1)=2;
end



% HERE i DO NOT PUT ANY CATOFF OF MANUAL REACTION TIMES
% TAKE ALL THE TRIALS, THE VERY SHORT AND THE VERY LONG.
% YOU CAN MAKE THE CUT OFF ANYTIME LATER
good_mrt= [Trials(good_trials).MRT]';
good_array=[Trials(good_trials).Array]';
good_trials_patterns=[Trials(good_trials).pattern]';

Good_tri=[good_trials good_mrt good_array good_trials_patterns];  %%%%%%%
%     Good_tri=[size(good_trials) size(good_mrt) size(good_array) size(good_trials_patterns)]

Good_tri(:,end+1)=1;
if ~isempty(Wrong_tri)
    all_trials=[Good_tri ;Wrong_tri];
else
        all_trials=  Good_tri;
end
all_sort=sortrows(all_trials,1);



%% CREATE THE LIST OF TRIALS FOR EACH CONDITION
for t=1:size(all_sort,1)
    type_trials(t,1)=Trials(all_sort(t,1)).Array;
    type_trials(t,2)=Trials(all_sort(t,1)).obj(1).pattern;
    type_trials(t,3)=Trials(all_sort(t,1)).obj(1).red;
    type_trials(t,4)=Trials(all_sort(t,1)).obj(1).green;
    type_trials(t,5)=Trials(all_sort(t,1)).obj(1).blue;
end
% % % %% check events code
% % % for xx=1:size(Good_tri,1)
% % %    
% % %     goodmatrr(xx,1)=Good_tri(xx,1);
% % %     moo =  [Trials(Good_tri(xx)).Events.Code];
% % %     event_no = find(moo(1,:) == 1012);
% % %     if ~isempty(event_no)
% % %         goodmatrr(xx,2)=event_no;
% % %     else
% % %         goodmatrr(xx,2)=NaN;
% % %     end
% % % 
% % % end
% % % 
% % % for xx=1:size(Wrong_tri,1)
% % %    
% % %     wrongmatrr(xx,1)=Wrong_tri(xx,1);
% % %     moo =  [Trials(Wrong_tri(xx)).Events.Code];
% % %     event_no = find(moo(1,:) == 1013);
% % %     if ~isempty(event_no)
% % %         wrongmatrr(xx,2)=event_no;
% % %     else
% % %         wrongmatrr(xx,2)=NaN;
% % %     end
% % % 
% % % end


%% divide two hands
% % % 
% % % numbcod=length(dat(2).Units);
% % % %% elimina trials
% % % %  all_sort = all_sort(20:end,:)
% % % 
% % % %% correct
% % % 
% % % 
% % % list1= all_sort(find(all_sort(:,5)==1),:);
% % % 
% % % %% Wrong
% % % list2= all_sort(find(all_sort(:,5)==2),:);
% % % 
% % % %% DIVIDE LEFT VS RIGHT HAND
% % % %% CORRECT and wrong insieme
% % % LH= list1(find(list1(:,3)==1),:);
% % % RH= list1(find(list1(:,3)==0),:);
% % % % %% WRONG
% % % 
% % % WRH= list2(find(list2(:,3)==1),:); % dovrebbe rispondee con la sinistra ma risponde con la destra
% % % WLH= list2(find(list2(:,3)==0),:);% dovrebbe rispondee con la destra ma risponde con la sinistra
% % % %% prendi solo mano destra
% % % %%sortrows(all_trials,1);
% % % all_RH=[];
% % % all_RH=[RH ;WRH];
% % % all_RH=sortrows(all_RH,1);
% % % 
% % % %%sortrows(all_trials,1);
% % % all_LH=[];
% % % all_LH=[LH ;WLH];
% % % all_LH=sortrows(all_LH,1);





%% SIX Some conversions -------
%  fn=[filenm '_spk']
% len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len));
% cd(PathName1);

clear Infos
Infos = NaN(size(Trials,1),20);
for i=1:size(dat,2)
%     i
    clear Temp;
     Temp = [[double(dat(i).Events.Code)];[double(dat(i).Events.Time)]]';
    
    if ~isempty(Temp(find((Temp(:,1)==800))))   & ~isempty(Temp(find((Temp(:,1)==1003)))) & ...
            ~isempty(Temp(find((Temp(:,1)==1007))))  & ~isempty(Temp(find((Temp(:,1)==1014)))) & ...
            ~isempty(Temp(find((Temp(:,1)==801))))   & ...
            (~isempty(Temp(find((Temp(:,1)==1012)))) | ~isempty(Temp(find((Temp(:,1)==1013)))))
%         
        
        Infos(i,1)  = Temp(find(Temp(:,1)== 800),2); %Start
        Infos(i,2)  = Temp(find(Temp(:,1)==1095),2); %Start Spike Recording
        Infos(i,3)  = Temp(find(Temp(:,1)==1003),2); %Monkey starts to fixate
        if size(Temp(find(Temp(:,1)==1007),2),1)>1
            bugvar=Temp(find(Temp(:,1)==1007),2);
            Infos(i,4)  = bugvar(1);
        else
            Infos(i,4)  = Temp(find(Temp(:,1)==1007),2); %Target
        end
        Infos(i,5)  = type_trials(i,2); %Target shape
        Infos(i,9)  = type_trials(i,1); %Hand (0-Right 1-Left)
        Infos(i,10) = all_sort(i,5);    %Correct = 1; Wrong = 2
%         Infos(i,10) = CORR;    %Correct = 1; Wrong = 2
        Infos(i,11) = Temp(find((Temp(:,1)==1012) | (Temp(:,1)==1013)),2);  %Time of correct or wrong
        Infos(i,12) = Temp(find(Temp(:,1)==1014),2); %End of trial
        Infos(i,13) = Temp(find(Temp(:,1)== 801),2); %End of recording
        Infos(i,14) = Infos(i,11)-Infos(i,4);    %RT
        Infos(i,15) = Infos(i,13)-Infos(i,1);    %Trial duration
        Infos(i,16) = 0;    %Category
        
           clear TEMP;
        if ~isempty(Temp(find(Temp(:,1)==1050),2)) & ~isempty(Temp(find(Temp(:,1)==1051),2))
            TEMP = Temp(find(Temp(:,1)==1051),2)-Temp(find(Temp(:,1)==1050),2);
            Infos(i,17) =  TEMP(find(TEMP>0));
        end
        
    else 
        
        
        LICK{i,1} = [];  
     end
end





disp('!!! END OF CODE CER_analysis_01_ADAPTEDL_n !!!');
disp('!!! Data ready for further analysis !!!');


CODE1_Pass = 1;




% end