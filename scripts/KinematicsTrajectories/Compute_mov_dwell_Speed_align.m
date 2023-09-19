function  [IntTrj_mov,IntTrj_dwell,get_mov_Speed,get_dwell_Speed] = Compute_mov_dwell_Speed_align(temp_Traj,params)
IntTrj_dwell= params.IntTrj_dwell;
InTrj_mov = params.InTrj_mov;
TapsBins = params.TapsBins;
seq = params.seq;
tap_thresh = params.tap_tresh;
mov = params.mov; 
mov_i = params.mov_i;
dwell = params.dwell; 
dwell_i = params.dwell_i;
cBin = params.cBin;
%% Properties

dur_mov = length(mov : mov_i)+5;%5 3
dur_dwell = length(dwell : dwell_i);

bw = floor(dur_mov/2); 
fw = dur_mov - bw;

if(cBin == 0 && params.seq ~= 1)

Start = TapsBins(1,seq) - bw+2 ;% - bw ; %2 para la figura
% Start = TapsBins(1,seq) +1 ;% - bw ; %2 para la figura

Tap_i_1 = TapsBins(1,seq) + fw;
dwell_i_1 = Tap_i_1 + dur_dwell-3; %- 1 -3
else
    IntTrj_mov = nan(1,InTrj_mov);
    IntTrj_dwell = nan(1,IntTrj_dwell);
    get_mov_Speed = nan;
    get_dwell_Speed = nan;
    return;
end

Traj_trace_all = temp_Traj(1:dwell_i_1);
Traj_trace = temp_Traj( Start :Tap_i_1);%+6
% Traj_trace = temp_Traj( Start :Tap_i_1 + 3);%+6

dwell_trace = temp_Traj( Tap_i_1 :dwell_i_1)';


IntTrj_mov = InterpspeedTraceTrajectory(Traj_trace',InTrj_mov);
IntTrj_dwell = InterpspeedTraceTrajectory(dwell_trace,IntTrj_dwell);

if(seq ~= 5)

get_mov_Speed = mean(Traj_trace); 
get_dwell_Speed = mean(dwell_trace);
else
get_mov_Speed = mean(Traj_trace); 
get_dwell_Speed = nan;
end

if(params.plot)
f2 = params.fig2;
set(0,'CurrentFigure',f2)
subplot(2,3,seq)
hold on;
plot(Traj_trace_all);
ymax = get(gca,'YLim');
Tapp = [TapsBins(1,seq),TapsBins(1,seq)];
movs = [Start Tap_i_1];
dwells = [Tap_i_1 dwell_i_1];
line([movs; movs],[ymax(1) ymax(2)],'LineWidth',1.5,'Color',[1 1 1])
plot([Tapp; Tapp],[ymax(1) ymax(2)],'c--')

if(seq ~= 5)
   plot([dwells; dwells],[ymax(1) ymax(2)],'r--')
 end
end
% close


