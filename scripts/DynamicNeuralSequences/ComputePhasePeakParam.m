%circplot is a function to plot circular statistics
%Input: temdata = data, maxScale = interval, modal=modality 1 for aud 2 for
%visual
function [outMeanR] = ComputePhasePeakParam(tPeakOnset,Param2, TeoDur,condi,SO,ploing,PlotOrder,tit,TapTime)

%figure setup
if ploing
    colormap('jet');
    cmap = colormap;
    MIL(1,:) = 'A450 ';
    MIL(2,:) = 'A850 ';
    MIL(3,:) = 'V450 ';
    MIL(4,:) = 'V850 ';
end
tPeakOnset(tPeakOnset == 99999) = [];
Param2(Param2 == 99999) = [];
PeakOnset = tPeakOnset - TapTime;

temdata = (PeakOnset.*360)./TeoDur; %Results in deegrees
temdata(temdata> 359.999999) = 360;
valsRad = deg2rad(temdata); %Results in radians

Nbins = 16;
stepedges = 2*pi/((Nbins));
edgesPolar = 0:stepedges:2*pi;
edgesPolar2 = edgesPolar;
edgesdegrees = rad2deg(edgesPolar);
N = histc(valsRad,edgesPolar);
N(16) = N(16)+N(17);
edgesPolar(end) = [];
N(end) = [];
hiscoun = N';
 NumberofCellsPerPhase = hiscoun;

% [N,edges] = histcounts(valsRad,Nbins);
% [NumberofCellsPerTimeBin,edges] = histcounts(PeakOnset,Nbins);

steptime = TeoDur/((Nbins));
Bintimes = 0:steptime:TeoDur;

% Bintimes = edges(1:end-1)';

maxhi = max(hiscoun);
temcoun = hiscoun/maxhi;
temcoun = temcoun';

if ploing
    scrsz = get(0,'ScreenSize');
  %  tit = ['Mean R Peak Phase'  ];
    f1 = figure(condi+PlotOrder);
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'white','Name', tit);
    whitebg([1 1 1])
    
    subplot(2,2,SO);
end

start = 1;
for ii = 1:16
    
    dat = NumberofCellsPerPhase(ii);
    ii
    r_anglet = circ_mean(valsRad(start:start+dat-1));
    resultantt = circ_r(valsRad(start:start+dat-1));
    numAct = numel(valsRad(start:start+dat-1));
    OutParam2 = nanmean(Param2(start:start+dat-1));
    outMeanR(ii,1) = ii;
    outMeanR(ii,2) = edgesPolar(ii);
    outMeanR(ii,3) = edgesdegrees(ii);
    outMeanR(ii,4) = Bintimes(ii);
    outMeanR(ii,5) = r_anglet;
    outMeanR(ii,6) = resultantt;
    outMeanR(ii,7) = NumberofCellsPerPhase(ii);
    outMeanR(ii,8) = OutParam2;
    
    theta = [edgesPolar(ii) edgesPolar(ii)];    
  %  theta = [r_anglet r_anglet];    
    rho = [0 OutParam2];
    
    if ploing
        polarplot(theta,rho,'LineWidth', 2, 'Color', [cmap(ii*3,:)]);
        hold on
    end
    start = start+dat;
end

if ploing
    title([MIL(condi,:) num2str(SO)],'FontSize',8);
    hold off
end

