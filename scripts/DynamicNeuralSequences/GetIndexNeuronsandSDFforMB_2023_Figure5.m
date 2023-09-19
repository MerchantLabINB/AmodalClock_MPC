Ttask = 4;
%load up last OutputStructFiltered Poission4Task;
temPoisson = Poission4Task;
DurModCellGroup = 12; %%all sign cells

%%%%Getting index of neurons per condition aout of 4
for task = 1:Ttask
    Neuidx = OutputStructFiltered{1,1}.IdxDurModSOET.DurModSOETIDs{task}(:,DurModCellGroup);
    Neuidx(Neuidx == 0) = [];
    TTneurons = numel(Neuidx);
    IndexNeurons4Condi{task} = Neuidx;
end

   
%%%getting the activiation period from poisson for the indexed neurons    
    for task = 1:Ttask
        Neuidx = IndexNeurons4Condi{task};
        TTneurons = numel(Neuidx);
        Flag = 1;
        
        for cellidx = 1:TTneurons
            
            cellidxOK = Neuidx(cellidx);
            clear SOActivations
            cellidxOK
            
            [SOActivations] = Get_DischargefromPoisson2021SO(temPoisson{cellidxOK,task},task);
            %TActivations{cellidx,task} = SOActivations;
            TActivations{cellidxOK,task} = SOActivations;
            if (SOActivations{1,1}.numac == 1 & Flag)
                TapTimes(task,1:5) =  SOActivations{1,1}.taptimes(:);
                Flag = 0;
            end
            
        end
    end
    
    
    TSDFs = [];
    BehTimesSDF = [];
        
    Neuidx = 1:1:numel(OutputStructFiltered{1,1}.Neurons);
    TTneurons = numel(Neuidx);
    
    for cellidx = 1:TTneurons
        
        cellidxOK = Neuidx(cellidx);
        cellidxOK
        
        [temSDFs,timessdf] = ComputeSDFs(Pop(cellidx,:));
        TSDFs = [TSDFs;temSDFs];        
        
    end    
    BehTimesSDF = [timessdf];
    
    
    save('DataMovingBumpAnalysis.mat','IndexNeurons4Condi','TActivations','TapTimes','TSDFs','BehTimesSDF');