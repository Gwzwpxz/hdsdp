addpath(genpath('.'));
cd src/potlp/potlp/mex/

try
    potlp_install;
catch
    fprintf("Installation failed :( \n");
end % End try

cd ../../../../
