function f = Simulation_Moving_Bumps()
% Figure 7. Neural sequence simulations and their resulting neural trajectories.

clear
f = figure('color', 'black','units','normalized','outerposition',[0 0 1 1]);
whitebg('k')

set(f, 'InvertHardcopy', 'off')
set(f, 'color', [0 0 0]);

% Modality  refers to absolute or relative
% Reset refers to resetting or not the neural sequences
% Sg_prop refers to 

% Panel A
Modality = 'Abs'; Reset = 'off'; Sg_prop = 'off'; Am_prop = 'off'; Motor = 'off'; Tuning = 'none';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 1)
% Panel B
Modality = 'Abs'; Reset = 'on'; Sg_prop = 'off'; Am_prop = 'off'; Motor = 'off'; Tuning = 'none';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 2)
% Panel C
Modality = 'Rel'; Reset = 'on'; Sg_prop = 'on'; Am_prop = 'off'; Motor = 'off'; Tuning = 'none';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 3)
% Panel D
Modality = 'Rel'; Reset = 'on'; Sg_prop = 'off'; Am_prop = 'on'; Motor = 'off'; Tuning = 'none';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 4)
% Panel E
Modality = 'Rel'; Reset = 'on'; Sg_prop = 'on'; Am_prop = 'off'; Motor = 'on'; Tuning = 'none';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 5)
% Panel F
Modality = 'Rel'; Reset = 'on'; Sg_prop = 'on'; Am_prop = 'off'; Motor = 'on'; Tuning = 'single';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 6)
% Panel G
Modality = 'Rel'; Reset = 'on'; Sg_prop = 'on'; Am_prop = 'off'; Motor = 'on'; Tuning = 'double';
plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, 7)


function plot_Sim(Modality, Reset, Sg_prop, Am_prop, Motor, Tuning, num_pan)

    Dur = [450 850]; % Durations
    nNeu = 400; % Number of neurons
    nDur = 2;   % Number of durations
    nSO = 2;    % Number f serial orders
    nB = 100;   % Number of bins in serial order
    sgp = 0.3;  % activation period scaling factor
    fcN = 50;   % number of shared neurons
    frc = 50;
    J = fcN+1:nNeu-fcN;
    
    for cDur = 1:2
        
        t{cDur} = linspace(0, nSO*Dur(cDur), nSO * nB + 1)';
        
        switch Modality
            case 'Abs'
                T = linspace(0, Dur(end), nNeu);
            case 'Rel'
                T = linspace(0, Dur(cDur), nNeu);    
        end
        
        
        switch Sg_prop
            case 'on'
                sg = sgp * Dur(cDur) * ones(1, nNeu);
            case 'off'
                sg = sgp * 650  * ones(1, nNeu);            
        end
    
        switch Motor
            case 'on'
                T(1:fcN) = linspace(0, Dur(cDur)/frc, fcN);
                T(nNeu-fcN+1:nNeu) = linspace((frc-1)*Dur(cDur)/frc, Dur(cDur), fcN);
                T(fcN+1: nNeu-fcN) = linspace(Dur(cDur)/frc + sgp * Dur(cDur), (frc-1)*Dur(cDur)/frc - sgp * Dur(cDur), nNeu - 2*fcN);
            case 'off'
                
        end
        switch Tuning
            case 'single'
                if cDur == 1
                    T(J(1:2:end)) = inf;
                end
            case 'double'
                if cDur == 1
                    T(J(1:2:end)) = inf;
                else
                    T(J(2:2:end)) = inf;
                end
            otherwise
                
        end
        switch Am_prop
            case 'on'
                A = Dur(cDur) / Dur(2);
            case 'off'
                A = 1;
        end
        
        switch Reset
            case 'on'
                d = mod(t{cDur}-T, Dur(cDur));
                d = min(d, Dur(cDur)-d);
            case 'off'
                d = abs(t{cDur}-T);
        end
        
        % spike density calculation
        R{cDur} = A * exp(-(d ./ sg).^2); 
        R{cDur}(:, T>Dur(cDur)) = 0;
        subplot(7, 4, 4*(num_pan-1)+cDur)
        % plot of spike density
        pcolor(t{cDur}, 1: nNeu, R{cDur}')
        caxis([0 1])
        axis square
        axis([0 2*Dur(end) 1 nNeu])
        shading interp
        set(gca, 'xtick', Dur(cDur)*(0:2), 'YTick', [1 nNeu])
        box off
        if cDur == 1
            ylabel('Neuron')
        end
        colormap turbo
    end
    
    % calculation of neural trajectories
    RT = [R{1};R{2}];
    [~, PC] = pca(normalize(RT));
    PC = PC(:, 1:3)';
    PM = mean(PC(:, [1 nSO*nB+2]), 2);
    ID = repmat(1:nDur, nSO*nB+1, 1); ID = ID(:);
    subplot(7, 4, 4*(num_pan-1)+3)
    hold on
    % colorplot = [0.75, 0, 0; 0.5, 0, 0];
    colorplot = [0,0.7490,1; 0 0 1];

    for c1 = 1: cDur    
        plot3(PC(1, ID == c1), PC(2, ID == c1), PC(3, ID == c1), 'LineWidth',3-c1, 'color', colorplot(c1, :))
        plot3(PC(1, (0:nB:nSO*nB)+(c1-1)*(nSO*nB)+c1), PC(2, (0:nB:nSO*nB)+(c1-1)*(nSO*nB)+c1), PC(3, (0:nB:nSO*nB)+(c1-1)*(nSO*nB)+c1), 'o', 'MarkerFaceColor',colorplot(c1, :))
    end
    plot3(PM(1), PM(2), PM(3), 'sk', 'MarkerFaceColor',[0 0 0])
    view(3)
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    axis equal
    set(gca, 'XTick', 0, 'YTick', 0, 'ZTick', 0)
    grid on
    hold off    
    ka = Amplitud_Index(PC, PM, nSO, nB, Dur);
    kv = Velocity_index(PC, nSO, nB, Dur);
    AMSI = AMSI_index(ka, kv);
    subplot(7, 4, 4*(num_pan-1)+4)
    plot(1: nSO*nB, AMSI, 'color','r','linewidth', 2)
    set(gca, 'tickdir', 'out','XTick', 0:nB:nSO*nB, 'XTickLabel', {'0' '\pi' '2 \pi'}, 'ytick', -1:1)
    grid on
    axis([1 nSO*nB+1 -1 1])
    axis square
    ylabel('AMSI')
end

function ka = Amplitud_Index(PC, PM, nSO, nB, Dur)
    fSig = 5;
    Amplitud = vecnorm(PC - PM);
    Amplitud = reshape(Amplitud, [nSO*nB+1, 2])';
    DiffAmp = diff(Amplitud);
    MinAmp = 0;
    MaxAmp = max(Amplitud, [], 'all') * (1-Dur(1)/Dur(2));
    ka = logsig(fSig*(-0.5+(DiffAmp-MinAmp) / (MaxAmp-MinAmp)));
end

function kv = Velocity_index(PC, nSO, nB, Dur)
    tR = linspace(0, nSO, nSO*nB+1);
    fSig = 5;
    Velocidad(1, :) = vecnorm(diff(PC(:, 1: nSO*nB+1), [], 2) ./ diff(Dur(1) * tR));
    Velocidad(2, :) = vecnorm(diff(PC(:, (1: nSO*nB+1) + nSO*nB+1), [], 2) ./ diff(Dur(2)*tR));
    DiffVel = diff(Velocidad([2 1], :));
    MinVel = 0;
    MaxVel = max(Velocidad, [], 'all') * (1-Dur(1)/Dur(2));
    kv = logsig(fSig*(-0.5+(DiffVel-MinVel) / (MaxVel-MinVel)));
end

function AMSI = AMSI_index(ka, kv)
    AMSI = (ka(:, 1:end-1) - kv) ./ (ka(:, 1: end-1) + kv);
end

end