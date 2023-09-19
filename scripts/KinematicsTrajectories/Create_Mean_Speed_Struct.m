function  trajStruct2 = Create_Mean_Speed_Struct(trajStruct)
for TaskCond  = 1:4
mov_ini = [];
mov_end = [];
dwell_i = [];
dwell_end = [];
dwell_pause = [];
tresh_seq = [];


for i=1:size(trajStruct,1)
    mov_ini = [mov_ini; trajStruct(i, TaskCond).mov_ini];
    mov_end = [mov_end; trajStruct(i, TaskCond).mov_end];
    dwell_i = [dwell_i; trajStruct(i, TaskCond).dwell_i ];
    dwell_end=  [dwell_end; trajStruct(i, TaskCond).dwell_end];
%     dwell_pause = [dwell_pause; trajStruct(i, TaskCond).dwell_pause];
    tresh_seq = [tresh_seq; trajStruct(i, TaskCond).bin_thresh_seq];  
end

% dwell_pause(dwell_pause == 0 ) = nan;
trajStruct2(1, TaskCond).mov_ini        = repmat(floor(mean(mov_ini)),25,1);
trajStruct2(1, TaskCond).mov_end        = repmat(floor(mean(mov_end)),25,1);
trajStruct2(1, TaskCond).dwell_i        = repmat(floor(mean(dwell_i)),25,1);
trajStruct2(1, TaskCond).dwell_end      = repmat(floor(mean(dwell_end)),25,1);
% trajStruct2(1, TaskCond).dwell_pause      = repmat(floor(mean(dwell_end)),25,1);
trajStruct2(1, TaskCond).bin_thresh_seq = repmat(floor(mean(tresh_seq)),25,1);
trajStruct2(1, TaskCond).trialId_Task   = 1:25;

end