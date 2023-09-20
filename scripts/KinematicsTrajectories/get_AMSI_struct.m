function  [Amplitude,Speed ] = get_AMSI_struct(kinematics)

Trajectories = kinematics.trajectories;


SpeedTrajInter = cellfun(@InterpTrajectory, Trajectories, 'UniformOutput', 0);
SpeedAmsi= cellfun(@Trajectory_Speed, SpeedTrajInter, 'UniformOutput', 0);

Amp_arr = reshape(kinematics.AmplitudeTrialsAMSI,[2,2])';

for i = 1:4
    AA = SpeedAmsi(:, i);
Speed_arr{1,i} = cat(2,AA{:})';
end
Speed_rep = reshape(Speed_arr,[2,2])';

for mod = 1:2
for dur = 1:2
Amplitude(1,:,dur,mod,:) = permute(Amp_arr{mod, dur},[2,1]);
Speed(1,:,dur,mod,:) = permute(Speed_rep{mod, dur},[2,1]);
end
end

Speed = cat(2, Speed(:, 1, :, :, :), Speed);