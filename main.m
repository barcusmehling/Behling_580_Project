%% Calculate flight forces while varying num. of accels and forces
clc;close all;clear all;
disp('This script takes a while to run...')
%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%
% nforces = 24, naccs = 20 for paper results
% These values can be varied, though the script has not been extensively
% tested for many other possible combinations...
nforces = 24;
naccs = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading Flight Environment...')
load Env_580_reduced.mat; % Simulated Accel. Flight Env - 135 x 135 x 401 PSD matrix (naccels x naccels x n freq. lines)
Y2 = Y2_new; % rename variable
disp('Loading Transfer Function Matrix...')
load FRF_580_reduced.mat; % Simulated Transfer / Frequency Response Function Matrix - 135 x 62 x 401 PSD matrix (naccels x nflightforces x n freq. lines)
H = H_new;

rng(42); % set random seed to get same results as poster
df = fs(2)-fs(1); % frequency spacing (fs imported from previous .mats)
rms_inds = 3:401; % calculate root-mean-square error ignoring first two freq. lines

errs_cval = zeros(naccs-1,nforces); % leave-one-out cross validation errors for each case
errs_global = zeros(naccs-1,nforces); % "global" error on uncontrolled accelerometers for each case

force_vec = randperm(62); % randomize order of model forces to be added one at a time

clc;

for zz = 2:naccs % loop through num. accel sets (no leave-one-out error for 1 accel case so start at 2 accel case)
    disp(['Running Accel Set ' num2str(zz) ' of ' num2str(naccs) '...'])
    inds = 1:zz; % include first "zz+2" accels 
    for ii = 1:nforces % increment through candidate forces
        force_inds = force_vec(1:ii); % retain first "ii" forces using randomized order
        err_cval = 0; % initialize cross validation error for this shaker-accel combo
        unc_ind = 1:size(Y2,1); % for "global error", calculate across all accels not included in training data
        unc_ind(inds) = []; % see previous line
        for jj = 1:length(inds) % calculate forces while leaving out one accel each time
            train_inds = inds; % training data indices
            err_ind = train_inds(jj); % index of left out accel to calc error (NOTE: assigning this value then modifying in next line would mess things up in python I think but works in MATLAB)
            train_inds(jj) = []; % remove left out accel from training data indices
            Sxxest = zeros(size(fs)); % estimated PSD matrix at left-out accel
            for kk = 1:length(fs) % calculate forces at each frequency line
                Sff = pinv(H(train_inds,force_inds,kk))*Y2(train_inds,train_inds,kk)*pinv(H(train_inds,force_inds,kk))'; % estimate force PSD 
                Sxxest(kk) = H(err_ind,force_inds,kk)*Sff*H(err_ind,force_inds,kk)'; % estimate env at left out accel
            end
            err_cval = err_cval + sqrt(sum((0.5*mag2db(abs(Sxxest(rms_inds)'))...
                -0.5*mag2db(squeeze(abs(Y2(err_ind,err_ind,rms_inds))))).^2)*df); % calculate RMS dB error at left out accel
        end
        % calculate "global" error by calculating forces when using all
        % training data points
        train_inds = inds; % training data indices
        Sxxfull = zeros(length(unc_ind),length(unc_ind),length(fs)); % estimated PSD matrix at "global" inds (all test inds)
        for kk = 1:length(fs) % estimate env
            Sff = pinv(H(train_inds,force_inds,kk))*Y2(train_inds,train_inds,kk)*pinv(H(train_inds,force_inds,kk))'; % estimate force PSD 
            Sxxfull(:,:,kk) = H(unc_ind,force_inds,kk)*Sff*H(unc_ind,force_inds,kk)'; % estimate full environment
        end
        [~,eval]= e_asd(fs,get_psd(abs(Sxxfull)),fs,get_psd(abs(Y2(unc_ind,unc_ind,:))),10,2000); % calculate error globally
        errs_cval(zz-1,ii) = err_cval / length(inds); % average sum to get leave-one-out error
        errs_global(zz-1,ii) = eval; % save global error
    end
end

% Define a set of unique line styles & markers
styles = {'-o','-s','-d','-^'};

plot_inds = [2 5 8 11]; % accel cases to plot

figure; % plot double-descent
hold on;
for ii = 1:length(plot_inds)
    plot(1:nforces, errs_cval(plot_inds(ii),:), styles{ii}, 'LineWidth', 2, 'MarkerSize', 6);
end
legend({'3 Acc','6 Acc','9 Acc','12 Acc'},'Location','southeast');
grid on;
xlabel('nforces');
ylabel('Leave-One-Out Error');
hold off;
xlim([0 16])

% Note: errors in the following plots are in different units technically,
% though the overall trend is of interest so it is ok
figure; % plot leave-one-out error contourmap for nforcs vs. naccels
contourf(errs_cval)
colorbar;
xlabel('nforces')
ylabel('naccels')
title('Leave-One-Out Error')

figure; % plot global error
contourf(errs_global)
colorbar;
xlabel('nforces')
ylabel('naccels')
title('Global Error')

clc;
disp('DONE')
%% Functions used in script

function psd = get_psd(Sxx) % retain diagonal terms only from PSD matrix
    psd = zeros(size(Sxx,1),size(Sxx,3));
    for ii = 1:size(Sxx,1)
        for jj = 1:size(Sxx,3)
            psd(ii,jj) = abs(Sxx(ii,ii,jj));
        end
    end
end

function [e_ASD,e_val] = e_asd(target_fs,target_psds,result_fs,result_psds,fmin,fmax)
    % This function written by Chris Schumann, UW-Madison, 2019
    % Calculates RMS dB error between two sets of PSDs

    % sort vectors to have same dimensions
    size_target = size(target_psds);
    if size_target(2) < size_target(1)
        target_psds = target_psds';
        size_target = size(target_psds);
    end
    size_result = size(result_psds);
    if size_result(2) < size_result(1)
        result_psds = result_psds';
        size_result = size(result_psds);
    end
    n_psds = size_result(1);
    
    % determine which frequency vector to interpolate to and do so
    if size_result(2) > size_target(2)
        fs1 = result_fs;
        psds1 = result_psds;
        fs2 = target_fs;
        psds2 = target_psds;
    else
        fs1 = target_fs;
        psds1 = target_psds;
        fs2 = result_fs;
        psds2 = result_psds;
    end
    
    % interpolate
    for jj = 1:n_psds
        psds2_int(jj,:) = interp1(fs2,psds2(jj,:),fs1);
    end
    
    % frequency range of interest
    [~,fmin_ind] = min(abs(fs1 - fmin));
    [~,fmax_ind] = min(abs(fs1 - fmax));
    
    % error calculate
    G_error_sum = zeros(length(fmin_ind:fmax_ind),1);
    for kk = 1:n_psds
        db_psds1 = 0.5*mag2db(psds1(kk,fmin_ind:fmax_ind));
        db_psds2 = 0.5*mag2db(psds2_int(kk,fmin_ind:fmax_ind));
        G_error(:,kk) = (db_psds1 - db_psds2).^2;
        G_error_sum = G_error(:,kk) + G_error_sum;
    end
    
    % output error at frequency lines
    e_ASD = sqrt(G_error_sum/n_psds);
    
    % output error value
    e_val = sqrt(sum(e_ASD.^2)/(length(fmin_ind:fmax_ind)));

end

