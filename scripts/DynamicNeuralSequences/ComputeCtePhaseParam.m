%circplot is a function to plot circular statistics
%Input: temdata = data, maxScale = interval, modal=modality 1 for aud 2 for
%visual
function [outMeanR] = ComputeCtePhaseParam(tPeakOnset,Param2, TeoDur,condi,SO,ploing,PlotOrder,tit)

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
PeakOnset = tPeakOnset - tPeakOnset(1);

temdata = (PeakOnset.*360)./TeoDur; %Results in deegrees
valsRad = deg2rad(temdata); %Results in radians
Nbins = 16;
[N,edges] = histcounts(valsRad,Nbins);
hiscoun = N';
stepedges = 2*pi/(numel(edges)-1);
edgesPolar = 0:stepedges:2*pi;
NumberofCellsPerPhase = hiscoun;
[NumberofCellsPerTimeBin,edges] = histcounts(PeakOnset,Nbins);



if ploing
    scrsz = get(0,'ScreenSize');
  %  tit = ['Mean R Peak Phase'  ];
    f1 = figure(condi+PlotOrder);
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'white','Name', tit);
    whitebg([1 1 1])
    
    subplot(2,2,SO);
end


    step = round(numel(valsRad)/16);
    temstep = numel(valsRad) - step*16;    
    laststep = step+temstep;
    
    if laststep  < 0
        step = step-1;
        temstep = numel(valsRad) - step*16;
        laststep = step+temstep;
    end
        
    dat = 1;
    
    for ii = 1:15
      
        r_anglet = circ_mean(valsRad(dat:(dat+step-1)));
        resultantt = circ_r(valsRad(dat:(dat+step-1)));
        numAct = numel(valsRad(dat:(dat+step-1)));
        outMeanR(ii,1) = ii;
        outMeanR(ii,2) = edgesPolar(ii);
        outMeanR(ii,3) = r_anglet;
        outMeanR(ii,4) = resultantt;
        outMeanR(ii,5) = numAct;
        theta = [r_anglet r_anglet];
        rho = [0 resultantt];
        
        if ploing
            polarplot(theta,rho,'LineWidth', 2, 'Color', [cmap(ii*3,:)]);
            hold on
        end
        
         dat = dat+step;%(ii-1) * step;
        
    end
    ii = 16;
    r_anglet = circ_mean(valsRad(dat:(dat+laststep-1)));
    resultantt = circ_r(valsRad(dat:(dat+laststep-1)));
    numAct = numel(valsRad(dat:(dat+laststep-1)));
    outMeanR(ii,1) = ii;
    outMeanR(ii,2) = edgesPolar(ii);
    outMeanR(ii,3) = r_anglet;
    outMeanR(ii,4) = resultantt;
    outMeanR(ii,5) = numAct;
    theta = [r_anglet r_anglet];
    rho = [0 resultantt];
    
    if ploing
        polarplot(theta,rho,'LineWidth', 2, 'Color', [cmap(ii*3,:)]);
        hold on
    end
    
 

if ploing
    title([MIL(condi,:) num2str(SO)],'FontSize',8);
    hold off
end

