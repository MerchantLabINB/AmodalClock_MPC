function CombineFigures(cross_correlograms,namefig)

% cross_correlograms =pwd;
% savepath = '.\';

% get contents of each folder
corrs = dir(fullfile(cross_correlograms, 'Panel*png'));

i = 1;
PanelA = imread(fullfile(corrs(i).folder, corrs(i).name));   
PanelB = imread(fullfile(corrs(i+1).folder, corrs(i+1).name));   
PanelC = imread(fullfile(corrs(i+2).folder, corrs(i+2).name));   
PanelD = imread(fullfile(corrs(i+3).folder, corrs(i+3).name));   
PanelE = imread(fullfile(corrs(i+4).folder, corrs(i+4).name));   
PanelFGH = imread(fullfile(corrs(i+5).folder, corrs(i+5).name));   

%Fill with Black Image
%For imgCorr2
[m,n,~]= size(PanelD);
suppPanel = zeros(m,(n*3 - n*2)/2,3);
suppPanel2 = zeros(m,(n*3 - n*2),3);

imgCorr1 = [PanelA PanelB PanelC];
imgCorr2 = [suppPanel PanelD PanelE suppPanel ];
imgCorr3 = [suppPanel2 PanelFGH suppPanel2];
corr_img = cat(1, imgCorr1, imgCorr2,imgCorr3);
f1 = figure;
f1.WindowState = 'maximized';
imshow(corr_img)
print(gcf,fullfile(cross_correlograms,namefig),'-dpdf','-r300','-fillpage')
imwrite(corr_img, fullfile(cross_correlograms,[namefig,'.tiff']));
close
%       imwrite(corr_img, fullfile(pwd, 'Figure2-PCA.tif'))
% % TODO: make this adjustment not necessary (or more automatic)
%     % these values cut off the margin around the figure
%     cutoff1 = round(size(corr_img, 2) * 0.09); 
%     cutoff2 = round(size(corr_img, 2) * 0.0741);
%     corr_img = corr_img(:,cutoff1:end - cutoff2,:);
% %     
% %     cutoff1 = round(size(PanelFGH, 2) * 0.105); 
% %     cutoff2 = round(size(PanelFGH, 2) * 0.08);
% %     PanelFGH = PanelFGH(:,cutoff1:end - cutoff2,:);
% 
% % we need to shrink the template figures so that we can append them to
%     % the cross_correlogram figure
%     % assuming ks neurons -> rows and lab_neurons -> columns in cross correlograms
% %     resize_factor_ks = (size(corr_img, 1) / length(PanelB)) / size(PanelC, 1); 
%       %-----Up Panel-----%
%       corr_img = cat(2, corr_img, PanelB,PanelC);
%       corr_img2 = cat(2, PanelD, PanelE);
%      
%       %-----Middle paneles-----%
%       binDifference= size(corr_img,2) - size(corr_img2,2);
%       addSamples = floor((size(corr_img,2) - size(corr_img2,2))/2);
%       addSamplesfw= binDifference - addSamples;
% %       [m,n,o]= size(corr_img,2)
%       addSpacebw =  zeros(size(corr_img,1), addSamples, size(corr_img,3));
%       addSpacefw =  zeros(size(corr_img,1), addSamplesfw, size(corr_img,3));
%       corr_img2 = cat(2, addSpacebw, corr_img2,addSpacefw);
% 
%       %-----Bottom Panel-----%
%      %Vertical
% %     binDifferenceV= size(PanelB,1) - size(PanelFGH,1);
% %     addSamples = floor((size(PanelB,1) - size(PanelFGH,1))/2);
% %     addSamplesfw= abs(binDifferenceV - addSamples);
% %     addSpacebw=  zeros(addSamples,size(PanelB,2), size(PanelB,3));
% %     addSpacefw =  zeros(addSamplesfw,size(PanelB,2), size(PanelB,3));
% %     PanelFGH = cat(1, addSpacebw, PanelFGH,addSpacefw);
% %-----Horizontal-----%
%     binDifference= size(corr_img,2) - size(PanelFGH,2);
%     addSamples = floor((size(corr_img,2) - size(PanelFGH,2))/2);
%     addSamplesfw= binDifference - addSamples;
%     addSpacebw=  zeros(size(corr_img,1), addSamples, size(corr_img,3));
%     addSpacefw =  zeros(size(corr_img,1), addSamplesfw, size(corr_img,3));
%     
%     corr_img3 = cat(2, addSpacebw, PanelFGH,addSpacefw);
% 
%     
%     
%     
%     
%       corr_img = cat(1, corr_img, corr_img2,corr_img3);
%       imwrite(corr_img, fullfile(pwd, 'Figure2-PCA.tif'));
% end