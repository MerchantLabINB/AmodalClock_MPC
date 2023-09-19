function SpeedParams = Get_trace_mean_Speed(trajectoriestoSpeed,trajStruct,lag)


trajectories = trajectoriestoSpeed;
Trajectory = trajectories;
%
% %Speed of Trajecties
SpeedTraj_temp = cellfun(@Trajectory_Speed, Trajectory, 'UniformOutput', 0);
% SpeedTraj = trajectoriestoSpeed;
if(lag==1)
    ControlBin = [0,0,0,0];
    idxNumlag = 6;%5,8
    SpeedTraj = cellfun(@(x) [nan(idxNumlag,1); x(1:end-idxNumlag)],SpeedTraj_temp,'un',0);
    %Para la correlación toma después del 9no bin.
    %     SpeedTraj = trajectoriestoSpeed;
else
    SpeedTraj = SpeedTraj_temp;
    ControlBin = [10,10,10,10];
end



if(trajStruct(1).avgSpeed)
    Speed_Modality = arrayfun(@(col) horzcat(SpeedTraj{:, col})', 1:size(SpeedTraj, 2), 'un', 0);
    Speed_Mean =cellfun(@(x) mean(x)',Speed_Modality,'un',0);
end

%---Article parameters figures-----%
% mean_movs_bin = [outputStruct.Movs];     %30    32    31    34
% mean_dwells_bin = [outputStruct.dwells]; %55   126    53   139
mean_movs_bin = repmat([30,32],1,2);
mean_dwells_bin = repmat([55,139],1,2);

%EDITED BY ABRAHAM 25/02/2023
%For article analysis
% mean_movs_bin = repmat([42,42],1,2);
% mean_dwells_bin = repmat([42,42],1,2);

% ControlBin = [10,10,10,10];
% ControlBin = [10,10,10,10];

%Trajectories Monkey 2
TapBins(1,:) = [1, 22:22:22*6] + ControlBin(2);%25
TapBins(2,:) = [1, 42:42:42*6] + ControlBin(3);%25
TapsBins = repmat(TapBins,2,1);

% %Trajectories M1M2 - Projection 2023
% TapsBins(1,:) = [44,82,121,159,197];
% TapsBins(2,:) = [1,(42:42:168)]+25;
% TapsBins(3,:) = [44,82,121,159,197];
% TapsBins(4,:) = [1,(42:42:168)]+25;




% params.plot = false;

for TaskCond = 1:4
    for Session = 1:size(trajStruct,1)

        IntTrj_mov = [];
        IntTrj_dwell = [];
        get_mov_Speed = [];
        get_dwell_Speed = [];
        lentrials = length((trajStruct(Session, TaskCond).trialId_Task));
        for trial = 1:lentrials
            params.plot = false;

            if(params.plot)
                f2  = figure(trial);
                whitebg('k')
                set(f2, 'InvertHardcopy', 'off')
                set(f2, 'color', [0 0 0]);
                params.fig2 = f2;
            end

            if(trajStruct(1).avgSpeed)
                temp_Traj = Speed_Mean{1,TaskCond};
            else
                temp_Traj = SpeedTraj{trial, TaskCond};
            end

            for seq = 1:5

                bin_thresh_seq = trajStruct(Session, TaskCond).bin_thresh_seq(trial,seq);
                params.TapsBins = TapsBins(TaskCond,:);
                params.InTrj_mov = mean_movs_bin(TaskCond);
                params.IntTrj_dwell = mean_dwells_bin(TaskCond);

                params.seq = seq;
                params.tap_tresh = bin_thresh_seq;
                params.mov = trajStruct(Session, TaskCond).mov_ini(trial,seq);
                params.mov_i =trajStruct(Session, TaskCond).mov_end(trial,seq);
                params.cBin = ControlBin(1,TaskCond);
                if(seq~=5)

                    %                     params.dwell = trajStruct(Session, TaskCond).mov_end(trial,seq);
                    %                     params.dwell_i = trajStruct(Session, TaskCond).mov_ini(trial,seq+1)-1;

                    params.dwell = params.mov_i;%trajStruct(Session, TaskCond).dwell_i(trial,seq);
                    params.dwell_i = trajStruct(Session, TaskCond).dwell_end(trial,seq);
                end
                %
                % [IntTrj_mov(seq,:,trial),IntTrj_dwell(seq,:,trial),get_mov_Speed(trial,seq),get_dwell_Speed(trial,seq)]=Compute_mov_dwell_Speed(temp_Traj,params);
                try
                    [IntTrj_mov(seq,:,trial),IntTrj_dwell(seq,:,trial),get_mov_Speed(trial,seq),get_dwell_Speed(trial,seq)]= Compute_mov_dwell_Speed_align(temp_Traj,params);
                catch
                end
            end
            %             close
        end

        %         %PAra generar el arreglo
        %         aux_mov = reshape(permute(IntTrj_mov,[3,2,1]),25,210);
        %         aux_dwell = reshape(permute(IntTrj_dwell,[3,2,1]),25,210);
        %         array_mov{TaskCond} = aux_mov(:,1:168);
        %         array_dwell{TaskCond} = aux_dwell(:,1:168);

        aux_dwell_mean = get_dwell_Speed(:,1:4);
        SpeedParams.mean_trace_mov{TaskCond}(Session,:) = mean(mean(IntTrj_mov,1,'omitnan'),3,'omitnan');
        SpeedParams.mean_trace_dwell{TaskCond}(Session,:) = mean(mean(IntTrj_dwell,1,'omitnan'),3,'omitnan');
        SpeedParams.mean_Mov(Session,TaskCond) = mean(get_mov_Speed(:),'omitnan');
        SpeedParams.mean_dwell(Session,TaskCond) = mean(aux_dwell_mean(:),'omitnan');
        SpeedParams.mov_SO{Session,TaskCond} = get_mov_Speed;
        SpeedParams.dwell_SO{Session,TaskCond} = aux_dwell_mean;

        %Modifed by Abraham Betancourt Vera
        %15 marzo 2023 for Canonical Correlation Analysis
        speed_dwell_segment = reshape(permute(IntTrj_dwell(1:4,:,:),[3,2,1]),25,[]);
        speed_mov_segment = reshape(permute(IntTrj_mov(1:4,:,:),[3,2,1]),25,[]);

        clear params
    end
    SpeedParams.dwell_traces{TaskCond}= speed_dwell_segment;
    SpeedParams.movs_traces{TaskCond} = speed_mov_segment;

end

