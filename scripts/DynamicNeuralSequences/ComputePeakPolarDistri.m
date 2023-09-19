%circplot is a function to plot circular statistics
%Input: temdata = data, maxScale = interval, modal=modality 1 for aud 2 for
%visual
function [OutCircStatistics] = ComputePeakPolarDistri(PeakOnset,TeoDur,BinPeak,halfPSYc,condi,ploing,PlotOrder)

%figure setup
if ploing
    ColorDose = [2 48];
    colormap('jet');
    cmap = colormap;
    ColorStep = [cmap(ColorDose(2),:)];
    MIL(1,:) = 'A450';
    MIL(2,:) = 'A850';
    MIL(3,:) = 'V450';
    MIL(4,:) = 'V850';
end


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

steptime = TeoDur/((Nbins));
Bintimes = 0:steptime:TeoDur;

maxhi = max(hiscoun);
temcoun = hiscoun/maxhi;
temcoun = temcoun';

[pvalRalegighAsyn] = circ_rtest(edgesPolar',temcoun'); %compute Raleigh Test
r_angle = circ_mean(edgesPolar',temcoun);
resultant = circ_r(edgesPolar',temcoun);
[~,circSTD] = circ_std(edgesPolar',temcoun);



OutCircStatistics.pvalRalegighAsyn = pvalRalegighAsyn;
OutCircStatistics.zRaleighAsyn = 0;
OutCircStatistics.r_angle = r_angle;
OutCircStatistics.resultant = resultant;
OutCircStatistics.circSTD = circSTD;


if ploing
    scrsz = get(0,'ScreenSize');
    tit = ['Peak Phase Histogram'  ];
    f1 = figure(condi+PlotOrder);
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'white','Name', tit);
    whitebg([1 1 1])
    
    subplot(2,2,condi);
    polarhistogram('BinEdges',edgesPolar2,'BinCounts',temcoun,...
        'FaceAlpha', 0.3, 'FaceColor',ColorStep, 'EdgeColor', 'none'); hold on
    
    theta = [r_angle r_angle];
    rho = [0 resultant];
    polarplot(theta,rho,'LineWidth', 2, 'Color',ColorStep)
    
    
    temPeak = (BinPeak.*360)./TeoDur; %Results in deegrees
    PeakRad = deg2rad(temPeak); %Results
    theta2 = [PeakRad PeakRad];
    rho2 = [0 0.5];
    
    ColorStep = [cmap(15,:)];
    polarplot(theta2,rho2,'LineWidth', 2, 'Color',ColorStep)
    
    
    temPSE = (halfPSYc.*360)./TeoDur; %Results in deegrees
    PSERad = deg2rad(temPSE); %Results
    theta3 = [PSERad PSERad];
    ColorStep = [cmap(35,:)];
    polarplot(theta3,rho2,'LineWidth', 2, 'Color',ColorStep);
    
    title([MIL(condi,:)],'FontSize',8);
    pax = gca;
    pax.FontSize = 12;
    pax.LineWidth = 2;
    angles = 0:45:360;
    pax.ThetaTick = angles;
    labels = {'0','','90','','180', '', '-90', ''};
    pax.ThetaTickLabel = labels;
    rlim([0 1]);  %rlim([0 150])
    rticks([0 0.25 0.5 0.75 1]);
    rticklabels({'0' '' '0.5' '' '1' })
    hold off
    %prepareForPrinting
end


if ploing
    scrsz = get(0,'ScreenSize');
    tit = ['Peak Phase Number Cells'  ];
    f1 = figure(condi+PlotOrder+5);
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'white','Name', tit);
    whitebg([1 1 1])
    
    subplot(2,2,SO);
end


start = 1;
for ii = 1:16
    
    dat = NumberofCellsPerPhase(ii);
 
    r_anglet = circ_mean(valsRad(start:start+dat-1));  
    numAct = numel(valsRad(start:start+dat-1));

    theta = [edgesPolar(ii) edgesPolar(ii)];   
  %  theta = [r_anglet r_anglet];
    rho = [0 numAct];
    
    if ploing
        polarplot(theta,rho,'LineWidth', 2, 'Color', [cmap(ii*3,:)]);
        hold on
    end
    start = start+dat;
end

if ploing
    title([MIL(condi,:)],'FontSize',8);
    hold off
end






