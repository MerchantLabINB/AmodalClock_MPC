function ax = AMSI(kinematics,f1,ax)
%Amplitud Modulation Temporal Scaling
%Create strcut for AMSI computation
[Amplitude,Speed ] = get_AMSI_struct(kinematics);

%ka - auditory
fSig = 5;
Dur = [450,850];
mAmp = reshape(mean(Amplitude,5),[168,2,2]);
Aamp  =mAmp(:,:,1)'; %2x168
DiffAmp = diff(Aamp);

MinAmp = 0;
MaxAmp = max(Amplitude, [], 'all') * (1-Dur(1)/Dur(2));
ka_a = logsig(fSig*(-0.5+(DiffAmp-MinAmp) / (MaxAmp-MinAmp)));

% Ka - Visual
fSig = 5;
Dur = [450,850];
mAmp = reshape(mean(Amplitude,5),[168,2,2]);
Aamp  =mAmp(:,:,2)'; %2x168
DiffAmp = diff(Aamp);

MinAmp = 0;
MaxAmp = max(Amplitude, [], 'all') * (1-Dur(1)/Dur(2));
ka_v = logsig(fSig*(-0.5+(DiffAmp-MinAmp) / (MaxAmp-MinAmp)));

%kv - auditory
mVel = reshape(mean(Speed,5),[168,2,2]);
DiffVel = diff(mVel(:,[2 1],1)');
MinVel = 0;
MaxVel = max(mVel, [], 'all') * (1-Dur(1)/Dur(2));
kv_A = logsig(fSig*(-0.5+(DiffVel-MinVel) / (MaxVel-MinVel)));

%kv - visual
mVel = reshape(mean(Speed,5),[168,2,2]);
DiffVel = diff(mVel(:,[2 1],2)');
MinVel = 0;
MaxVel = max(mVel, [], 'all') * (1-Dur(1)/Dur(2));
kv_V = logsig(fSig*(-0.5+(DiffVel-MinVel) / (MaxVel-MinVel)));

%AMSI
k_aud = (ka_a - kv_A) ./ (ka_a + kv_A);
k_visual = (ka_v - kv_V) ./ (ka_v + kv_V);

%% PLOT AMSI
set(0,'CurrentFigure',f1)

cols = 5;
rows = 8;

cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);
SP = 42:42:168;


ax(14) = subplot(rows,cols,[24,25]);
hold on;
plot(ka_v,'Color',cmap(2,:),'LineWidth',1.5)
plot(ka_a,'Color',cmap(4,:),'LineWidth',1.5)
ylim([0 1])
line([SP; SP],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--','LineWidth',1)

ax(15) =  subplot(rows,cols,[29,30]);
hold on;
plot(kv_V,'Color',cmap(2,:),'LineWidth',1.5)
plot(kv_A,'Color',cmap(4,:),'LineWidth',1.5)
% line([SP; SP],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--')
ylim([0 1])
line([SP; SP],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--','LineWidth',1)


ax(16)  = subplot(rows,cols,[34,35]);
hold on;
plot(k_aud,'Color',cmap(2,:),'LineWidth',1.5)
plot(k_visual,'Color',cmap(4,:),'LineWidth',1.5)
line([SP; SP],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--','LineWidth',1)


