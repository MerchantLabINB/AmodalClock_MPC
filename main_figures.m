run Startup

%Before you run this script, don't forget to download the
%pre-processed data and add them to the data folder.

%Figure 1: B,C,E
%% Modelo por orden serial
pathTofolder = '.\data';
pathToImages = '.\Figures';

if exist(pathTofolder,'dir')
    foldercontent = dir(pathTofolder);
    if numel(foldercontent) > 2
        sprintf('Warning: folder exists and is not empty:\t%s', 'data')


        %%
        %Figure 1 - Rhythmic tapping behavior.
        f1 = Behavior_analyses(pathTofolder);
        print(f1, '-dpdf', 'FigureS6_C.pdf','-bestfit');
        print(f1,fullfile(pathToImages,'Figure_1.png'),'-dpng','-r300');
        close(f1)
        %Figure 2 - Neural population trajectories.
        f2 = Trajectories_analyses(pathTofolder);
        print(f2,fullfile(pathToImages,'Figure_2.png'),'-dpng','-r300');
        close(f2)
        %Figure 3 - Kinematic of neural trajectories
        f3 = Kinematics_analysis(pathTofolder);
        print(f3,fullfile(pathToImages,'Figure_3.png'),'-dpng','-r300');
        close(f3)
        %Figure 4 - Properties of Neural Sequences.
        f4 = NeuralSequences_analyses(pathTofolder);
        print(f4,fullfile(pathToImages,'Figure_4.png'),'-dpng','-r300');
        close(f4)
        %Figure 5 - Dynamics of neural sequences.
        f5 = DNS_analyses(pathTofolder);
        print(f5,fullfile(pathToImages,'Figure_5.png'),'-dpng','-r300');
        close(f5)
        %Figure 6 - Generalizability of neural sequences.
        % ComputeDistanceDiagonalMatrices2023
        %%Figure 7 - . Neural sequence simulations.
        f7 = Simulation_Moving_Bumps();
        print(f7,fullfile(pathToImages,'Figure_7.png'),'-dpng','-r300');
       close(f7)

    else
        sprintf('Warning: you have to download data from Zenodo:\t%s', 'DOI:')
    end
end


