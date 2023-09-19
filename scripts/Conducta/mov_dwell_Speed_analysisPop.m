function  [movs_dwell,outputStruct,traj] = mov_dwell_Speed_analysisPop(Pop)
close all

%% Figure Properties
% f1 = figure(1);
% whitebg('k')
% f1 = gcf;
% set(f1, 'InvertHardcopy', 'off')
% set(f1, 'color', [0 0 0]);
% 
% params.fig2 = f2;


max_anorm_length=15; %40 
mov_threshold=20;%10  20
[m,n] = size(Pop);
params.plot = false;


bin_size = Pop(1).Times{1, 1}(2)-Pop(1).Times{1, 1}(1);


for TaskCond = 1:n
    arrayMov = [];
    arraydwells = [];
    dwells_array = [];
    movs_array = [];
   ArrayMovs =nan(25,5,14);
    Arraydwells =nan(25,4,14);
    for Session = 1:m

        idxTrial = Pop(Session, TaskCond).BehaveId;

        for trial = 1: numel(Pop(Session, TaskCond).Speed)
%             bin_size = Pop(TaskCond).Times{1, 1}(2)-Pop(TaskCond).Times{trial, 1}(1);
            
            disp(['Task ', num2str(TaskCond), ' Session ', num2str(Session),'..... Trial ', num2str(trial)])
            mov_dwell_array(trial).mov=[];
            mov_dwell_array(trial).dwell=[];
            mov_dwell_array(trial).mov_i=[];
            mov_dwell_array(trial).dwell_i=[];
            mov_dwell_array(trial).raw_vel= [];
            mov_dwell_array(trial).mov_dur = [];
            mov_dwell_array(trial).dwells = [];
            taps = Pop(Session, TaskCond).TimesTap(trial,2:end);
%             taps(end) = taps(end)-.150
            VideoSpeed = Pop(Session, TaskCond).Speed{trial, 1};
            VideoTime = Pop(Session, TaskCond).Times{trial,1};
            
            if(params.plot)
                whitebg('k')
                f1 = figure;
                set(f1, 'InvertHardcopy', 'off')
                set(f1, 'color', [0 0 0]);
            end
            
            movss = [];
            ddwells = [];
%             temp_Traj = SpeedTraj{trial, TaskCond};
            for seq = 1:5
                movs_aux = [];
                dwells_aux = [];
                mov_ii = [];

                numberofpulses = 1;
                ctrhes = 0;
                ctrhes2 = 0;
                flag = 0;
                cont=0;
%                   while((numberofpulses < 2) && (lenmovs <2))
                while(numberofpulses )
    
                    
                    if(flag)
                        mov_dwell_array(trial).mov(:)=[];
                        mov_dwell_array(trial).dwell(:)=[];
                        mov_dwell_array(trial).mov_i(:)=[];
                        mov_dwell_array(trial).dwell_i(:)=[];
                        mov_dwell_array(trial).raw_vel(:) = [];
                        movs_aux = [];
                        mov_ii = [];
                        dwells_aux = [];
                
                    end
                                 

                        thresh_seq2 = 0.01 + ctrhes;
                        thresh_seq = .20 +  ctrhes2 ;

                    
                    VideoTimeseries = VideoTime(VideoTime >= taps(seq)-thresh_seq & VideoTime <= taps(seq+1) + thresh_seq2 );%.100
                    
                    vel_cut = VideoSpeed(VideoTime >= taps(seq)-thresh_seq & VideoTime <= taps(seq+1) + thresh_seq2);%.100
                    mov_dwell_array(trial).videoTime = VideoTimeseries;
                    mov_dwell_array(trial).raw_vel=vel_cut;
                    vel_cut=vel_cut>mov_threshold;
                    mov_dwell_array(trial).vel_cut=vel_cut;
                    timems = mov_dwell_array(trial).videoTime*1000;
                    
                    c=0;
                    first=1;
                    vel_proc=[];
                    
                    
                    last_state=1;
                    for idx=1:length(vel_cut)
                        if(last_state==vel_cut(idx))
                            vel_proc(idx)=vel_cut(idx);
                            new_state=last_state;
                        else
                            lastIdx=idx+max_anorm_length;
                            if(lastIdx>length(vel_cut))
                                lastIdx=length(vel_cut);
                            end
                            
                            sum_states=sum(vel_cut(idx:lastIdx));
                            if(sum_states>5)%2
                                %stable state movement
                                vel_proc(idx)=1;
                                new_state=1;
                            else
                                %stable state dwell
                                vel_proc(idx)=0;
                                new_state=0;
                            end
                            
                            
                            
                        end
                        if(new_state==last_state)
                            c=c+1;
                        else
                            if(first<=2) %Remove first movement and first dwell
                                first=first+1;
                            else
                                if(last_state==1)
                                    movs_aux(end+1) = c; %#ok<AGROW>
                                    mov_ii(end+1) = idx;
                                    mov_dwell_array(trial).mov(end+1)=c;
                                    mov_dwell_array(trial).mov_i(end+1)=idx;
                                    
                                else
                                    dwells_aux(end+1)=c;
                                    mov_dwell_array(trial).dwell(end+1)=c;
                                    mov_dwell_array(trial).dwell_i(end+1)=idx;
                                    
                                end
                            end
                            c=0;
                        end
                        last_state=new_state;
                        %        mov_dwell_array(trial).tar_interval=fix(mean(intStims)*10);
                    end
                    
                    
                    
                    mov_dwell_array(trial).vector=vel_proc;
                    tempMov = mov_dwell_array(trial).mov;
                    lenMovs =  movs_aux;
                    %Find Pulses
                    numberofpulses = sum(diff(vel_proc*2 > 1 )  < 0);
                    
                    
                    if(numberofpulses > 3)      
                    ctrhes = ctrhes - .001;
                    flagpulses = 1;
                    elseif(numberofpulses == 1)
                    ctrhes = ctrhes + .001;
    
                    elseif(numberofpulses ~= length(dwells_aux)  && ~isempty(dwells_aux))
                       numberofpulses = 0;
                    elseif(numberofpulses == length(dwells_aux))
                    ctrhes = ctrhes + .0001;%-

%                     elseif(numberofpulses == length(dwells_aux) && flagpulses)
%                     ctrhes = ctrhes + .001;
                    else
                        
                     ctrhes = ctrhes + .01;
                    end
                    
                    if(cont > 1000)
                        ctrhes2 = ctrhes2 - .001;
                    end
                    flag = 1;
                    lenmovs = length(lenMovs);
                    cont = cont+1;
                end
                
%                 tap_threshold(trial,seq) = thresh_seq;

                %Show figure
                %Just for Article Panel
                if(params.plot)
                set(0,'CurrentFigure',f1)
%                 figure;
                subplot(2,3,seq)
%                 plot(timems,mov_dwell_array(trial).raw_vel)
%                 hold on; plot(timems,mov_dwell_array(trial).vector * 400)
                plot(mov_dwell_array(trial).raw_vel)
                hold on; plot(mov_dwell_array(trial).vector * 400)
                xlabel('Time (ms)')
                ylabel('Speed (a.u.)')
                title([' Mov ', num2str(seq)])
                end
     
           try
                              
           [movs2,ddwells2,mov_ini,mov_end] = Sort_dwell_movs(movs_aux,dwells_aux,mov_ii);
            
            movss  = [movss,movs2];
            ddwells = [ddwells, ddwells2];
            mov_end = (mov_ini + movs2);
            mov_ini_array(trial,seq) = round((mov_ini * (bin_size*1000))/20);
%             mov_end_array(trial,seq) =floor((mov_end * (bin_size*1000))/20);
            mov_end_array(trial,seq) =round(((mov_ini + movs2) * (bin_size*1000))/20) + 4;
            
            bin_thresh_seq(trial,seq) =  round((thresh_seq * 1000)/20);
            
            %% For speed Movement
            dwell_ini = mov_end;%3
            dwell_end = mov_end + ddwells2-1;%length(mov_end : mov_end + ddwells2-1);
%             dwell_pause(trial,seq) = round((ddwells_test* (bin_size*1000))/20);%length(mov_end : mov_end + ddwells2-1);

            dwell_bin_i(trial,seq) = round((dwell_ini * (bin_size*1000))/20);
            dwell_bin_end(trial,seq) = round((dwell_end * (bin_size*1000))/20)+4;
            

            catch e
               
           end
               
            end
            mov_dwell_array(trial).mov_dur =  movss;
            mov_dwell_array(trial).dwells = ddwells;
            
%             mov_dwell_array(idxTrial).mov_dur =  movss;
%             mov_dwell_array(idxTrial).dwells = ddwells; 
%         close
        end
                dwells= {mov_dwell_array.dwells};
                movs= {mov_dwell_array.mov_dur};   
                
                
                idxTrials = find(cellfun(@(x) length(x),movs,'un',1) == 5);
                movs = movs(idxTrials);
                dwells = dwells(idxTrials);
                
                                
                tempmov = cat(1,movs{:});
                temdwell = cat(1,dwells{:});
                
                
                Binary_movs = tempmov >= 20;
                true_movs = find(all(Binary_movs,2));

                tempmov = tempmov(true_movs,:);
                temdwell = temdwell(true_movs,:);

                %% Trajectories Params
                traj(Session,TaskCond).mov_ini = mov_ini_array(true_movs,:);
                traj(Session,TaskCond).mov_end = mov_end_array(true_movs,:);
                traj(Session,TaskCond).dwell_i = dwell_bin_i(true_movs,1:4);
                traj(Session,TaskCond).dwell_end = dwell_bin_end(true_movs,1:4);
                traj(Session,TaskCond).bin_thresh_seq = bin_thresh_seq(true_movs,:);
                traj(Session,TaskCond).trialId_Task= Pop(Session, TaskCond).BehaveId(true_movs)';  
                %%
                arrayMov = [arrayMov;tempmov];
                arraydwells = [arraydwells;temdwell];
                
%                 close all
               movs = tempmov * bin_size * 1000;
               dwells = temdwell * bin_size * 1000;
               
                ArrayMovs(idxTrial(true_movs),:,Session) = movs;
                Arraydwells(idxTrial(true_movs),1:4,Session) = dwells(:,1:4);
%                 auxdwellCell{TaskCond} = Arraydwells(:,1:4,Session);

%                movs = cell2mat(cellfun(@(x) (x()*bin_size)*1000,movs,'un',0));
%                 dwells = cell2mat(cellfun(@(x) (x()*bin_size)*1000,dwells,'un',0));
                
                
                dwells_array = [dwells_array;dwells(:)];
                movs_array =   [movs_array;movs(:)];
                
                AvgMovs(Session,TaskCond) =  mean(movs(:));
                Avgdwells(Session,TaskCond) = mean(dwells(:));
        clear mov_dwell_array 
%         close all
    end
    
        auxMovsCell{TaskCond} = ArrayMovs;
        auxdwellCell{TaskCond} = Arraydwells;

        outputStruct(TaskCond).Movs = floor(mean(arrayMov(:)));
        outputStruct(TaskCond).dwells = floor(mean(arraydwells(:)));
        
        movs_dwell.movs{TaskCond} =movs_array;
        movs_dwell.dwells{TaskCond} =dwells_array;
    
end

movs_dwell.avg_movs =AvgMovs;   
movs_dwell.avg_dwell =Avgdwells; 
movs_dwell.MovTrials = auxMovsCell;
movs_dwell.DwellTrials = auxdwellCell;
mean(dwells_array)
mean(movs_array)
disp('Revisa')









