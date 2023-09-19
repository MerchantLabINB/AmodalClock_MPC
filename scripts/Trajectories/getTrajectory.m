
function trajectory=getTrajectory(params,coeff,neuData)
    if(params.substractMean==1)
        neuData=bsxfun(@minus,neuData,mean(neuData,2));
    end
    if(~isempty(neuData))
        tic
        p=(coeff'*neuData);
%         p=(coeff(1:10,1:10)'*neuData(1:10,:));%Subspaces

        toc
%         col = p'; %Descomentar para no suavizar
        col=[]; %Comentar cuando no se suaviza

        %Denoise trajectory
        %-----Comentar para no suavizar-----%
        for i=1:size(p,1)
            col=horzcat(col,smooth(p(i,:),20,'loess'));
        end
        trajectory=col';
    else
        trajectory=[];
    end

end
