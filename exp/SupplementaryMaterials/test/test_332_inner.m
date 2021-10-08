% Inner test script for Synthetic Blind Deconvolution Minibatch (3.3.2)
% Experiment setup should be isolated in console file
% Variables already provided in outer environment
% 1. pfail 2. kappa 3. batchrange 4. steprange 5. data
U = data.U;
V = data.V;
b = data.b;
bestloss = data.bestloss;
tol = bestloss * 1.5;

[m, n] = size(U);

% Initialize result arrays
nsgdIterToOpt = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nproxlinIterToOpt = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;

% Do experiment with batchsize 1
batchsize = 1;

idx = 0;
for stepsize = steprange
    
    idx = idx + 1;
    tempSgdIter = ones(nTest, 1) * maxiter * m;
    tempProxlinIter = ones(nTest, 1) * maxiter * m;
    
    parfor i = 1:nTest
        
        init_z = randn(2 * n, 1);
        
        [sgdsol, sgdinfo] = proxsgdblind(U, V, b, sqrt(maxiter * m), 0, init_z, ...
            maxiter, tol, true, 0, stepsize, show_info);
        
        [proxlinsol, proxlininfo] = proxlinblind(U, V, b, sqrt(maxiter * m), 0, init_z, maxiter, tol, ...
            early_stop, 0, stepsize, show_info);
        
        if sgdinfo.status == "Optimal"
            tempSgdIter(i) = sgdinfo.niter;
        end % End if
        
        if proxlininfo.status == "Optimal"
            tempProxlinIter(i) = proxlininfo.niter;
        end % End if
        
    end % End parfor
    
    nsgdIterToOpt(1, idx) = mean(tempSgdIter);
    nproxlinIterToOpt(1, idx) = mean(tempProxlinIter);
    
end % End for

fprintf("- Batchsize 1 done \n");

for k = 1:nbatchtotest
    
    batchsize = batchrange(k);
    
    idx = 0;
    for stepsize = steprange
        
        idx = idx + 1;
        tempSgdIter = ones(nTest, 1) * maxiter * m;
        tempProxlinIter = ones(nTest, 1) * maxiter * m;
        
        parfor i = 1:nTest
            
            init_z = randn(2 * n, 1);
            
            [sgdsol, sgdinfo] = proxsgdblindbatch(U, V, b, sqrt(maxiter * m / batchsize), 0, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
            
            [proxlinsol, proxlininfo] = proxlinblindbatch(U, V, b, sqrt(maxiter * m / batchsize), 0, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
            
            if sgdinfo.status == "Optimal"
                tempSgdIter(i) = sgdinfo.niter;
            end % End if
            
            if proxlininfo.status == "Optimal"
                tempProxlinIter(i) = proxlininfo.niter;
            end % End if
            
        end % End parfor
        
        nsgdIterToOpt(k + 1, idx) = mean(tempSgdIter);
        nproxlinIterToOpt(k + 1, idx) = mean(tempProxlinIter);
        
    end % End for
    
    fprintf("- Batchsize " + batchrange(k) + " done \n");
    
end % End for

envname = "blind_kappa_" + kappa + "_pfail_" + pfail + "env.mat";
fprintf("Saving environment to " + envname + ". \n");
save(envname);

sgdSpeedup = min((nsgdIterToOpt'));
nsgdBatchOneIter = sgdSpeedup(1);
proxlinSpeedup = min((nproxlinIterToOpt'));
nproxlinBatchOneIter = proxlinSpeedup(1);
sgdSpeedup = max(sgdSpeedup) ./ sgdSpeedup;
proxlinSpeedup = max(proxlinSpeedup) ./ proxlinSpeedup;

fprintf("Experiments ended. Start plotting. \n");
fprintf("Type 1 figure: Speedup vs. Batchsize \n");

% Plot figure type 1: speedup
plot([1, batchrange], proxlinSpeedup, "-x",  "LineWidth", 3);
hold on;
plot([1, batchrange], sgdSpeedup, "-o", "LineWidth", 3);
legend(["SPL", "SGD"], "FontSize", 20);
% xlabel("batchsize m");
% ylabel("speedup");
set(gca, "Fontsize", 20);
savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_type_1.fig")
fprintf("Done. \n");

close all;

fprintf("Type 2 figure: Speedup vs. Stepsize \n");

% Plot figure type 2: robustness
for i = 1:nrobusttest
    
    idx = robustnessbatchidx(i) + 1;
    batchsize = robustnessbatch(i);
    
    fprintf("- Type 2 figure : Batchsize " + batchsize + " \n");
    
    robustsgdSpeedup = nsgdBatchOneIter ./ nsgdIterToOpt(idx, :);
    robustproxlinSpeedup = nproxlinBatchOneIter ./ nproxlinIterToOpt(idx, :);
    
    semilogx(steprange, robustproxlinSpeedup, "-x", "LineWidth", 3);
    hold on;
    semilogx(steprange, robustsgdSpeedup, "-o", "LineWidth", 3);
    
    legend(["SPL", "SGD"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_type_2_batch_" + batchsize + ".fig")
    
    %     xlabel("stepsize");
    %     ylabel("speedup");
    
    hold off;
    close all;
    
end % End

fprintf("Done. \n");

fprintf("Type 3 figure: Iteration to Opt vs. Stepsize");

% Plot figure type 3: Iteration to Opt
for i = 1:nrobusttest
    
    idx = robustnessbatchidx(i) + 1;
    batchsize = robustnessbatch(i);
    
    fprintf("- Type 3 figure : Batchsize " + batchsize + " \n");
    
    semilogx(steprange, nproxlinIterToOpt(idx, :), "-x", "LineWidth", 3);
    hold on;
    semilogx(steprange, nsgdIterToOpt(idx, :), "-o", "LineWidth", 3);
    
    legend(["SPL", "SGD"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_type_3_batch_" + batchsize + ".fig")
    
    %     xlabel("stepsize");
    %     ylabel("speedup");
    
    hold off;
    close all;
    
end % End

fprintf("Done. \n");

fprintf("**************************\n");
fprintf("*********  Done  *********\n");
fprintf("**************************\n");


