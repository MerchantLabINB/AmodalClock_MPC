%Run this file to add the directories to your matlab path
%Dowload the code and save AmodalClock_Merchant_Lab;
path = '..\AmodalClock_Merchant_Lab';%


%Crate Folder to save Figures
pathToImages = '.\Figures';
if ~exist(pathToImages, 'dir')
    mkdir(pathToImages)
    mkdir('.\data')

else
    sprintf('Warning: Folder  exist:\t%s', 'Figures')
end

addpath(genpath(fullfile(path,'data')))
addpath(genpath(fullfile(path,'scripts')))
addpath(genpath(fullfile(path,'Figures')))

%Before you run this file, you should download following
% repositories from MathWorks:
% boundedline.m
% platonic_solid.m
% Circular Statistics Toolbox
%TubPlot


