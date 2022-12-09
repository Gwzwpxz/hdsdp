clear; clc; close all;

cd /Users/gaowenzhi/Desktop/gwz/potred/matlab/test;
addpath ../src/;
addpath ../data/;

files = dir(fullfile("..", "data", "p_*"));
nfiles = length(files);
fnames = {files.name}';
maxmn = 1e+10;
minmn = 0;

% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.verbose = 1;
param.relFeasTol = 1e-04;
param.relOptTol = 1e-04;
param.QPWindow = 32;
param.RecordFreq = 100;
param.maxIter = 1000000;
param.maxRuizIter = 100;
param.maxPCIter = 0;
param.coefScal = 1;
param.curvature = 1000;
param.curvInterval = 500;
param.RScalFreq = 5;
param.PI_RestartMax = 10;
param.PI_RestartRate = 1.5;
param.maxTime = 600.0;
param.maxIter = 1000000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nfiles <= 1
    return;
end % End if

diary off;

log_path = "log_20221208_1";
mkdir(log_path);

fileID = fopen('1208_1.txt', 'w');

fieldnames = fieldnames(param);
nfields = length(fieldnames);
fprintf(fileID, "Parameter summary \n");
fprintf(fileID, "%-20s %5s \n", "Parameter", "Value");

for i = 1:nfields
    fprintf(fileID, "%-20s %f \n", fieldnames{i}, getfield(param, fieldnames{i}));
end % End for

fprintf(fileID, "| %30s |   pObj   |   dObj   | pInfeas | dInfeas | relGap  |  Time | Status \n",...
        "Instance");

for i = 1:nfiles
    
    n = fnames{i};
    dname = fullfile(log_path, n(3:end-4) + ".txt");
    
    if exist(dname, 'file')
        continue;
    end % End if
    
    diary(dname);
    diary on
    try 
        data = preprocess(fullfile("..", "data", fnames{i}));
        test_netlib(fullfile("..", "data", fnames{i}), param, maxmn, minmn, fileID);
    catch
        
    end % End try
    diary off;
end % End for
