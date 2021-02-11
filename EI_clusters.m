% Demo script for running spiking network simulations and analyses
% 
% by Luca Mazzucato 2020
%
% ----------------------------------------
% Please cite:
% L. Mazzucato, G. La Camera, A. Fontanini 
% Expectation-induced modulation of metastable activity underlies faster coding of sensory stimuli, 
% Nat. Neuro. 22, 787-796 (2019).
% ----------------------------------------
% This script runs simulations of LIF networks with excitatory (E) and
% inhibitory (I) spiking neurons, with E and I clusters  [from Wyrick et al. 2020]

% CUSTOM OPTIONS:
% 1) edit your custom-made stimuli inside auxi.create_params. 
%------------------------
ClustersOption='EI';%
%------------------------
% LOAD PARAMETERS
%------------------------
paramsfile='params.mat'; % file where all network parameters are saved
parameters(paramsfile);
savedir=fullfile('data'); mkdir(savedir); % setup directory for saving HMM data

%% RUN SIMULATION
ntrials=1; % number of trials
file_sim=fullfile(savedir,'results.mat');  % file where simulation results are saved
%---------------------------
% GENERATE SYNAPTIC WEIGHTS
%---------------------------
% J = N x N matrix of synaptic weights
% params = structure containing all network parameters
[J, params]=weightmatrix(paramsfile);
% [stimulus_save, params]=auxi.fun_stim(params); % STIMULUS

%% find threshold
ni=[params.ni_e params.ni_i]; % baseline firing rates in E and I populations from parameters.m
theta=auxi.fun_Fix_Threshold(ni,0,params);
params.theta_e=theta(1); 
params.theta_i=theta(2);
params.Theta=theta;
fprintf('  --- Fix_Threshold found:\n');
fprintf('      theta_e=%g; theta_i=%g \n',theta(1),theta(2));
fprintf('      Background rates: ni_e=%g, ni_i=%g\n',ni(1),ni(2));

%%
%------------------------
% SIMULATION
%------------------------
tic
firings=cell(1,ntrials); % IMPORTANT: cell array with all spike times in each trial. Each element is a 2-column array where the second column is the spike time, and the first column is the index of the neuron that fired
PlotData=cell(1,ntrials); % cell array with data for plotting
for iTrial=1:ntrials
    ParamsRun=params;
    ParamsRun.J=J;
    fprintf('--- Start SIM ...\n');
    [firings{iTrial}, PlotData{iTrial}]=simulation(ParamsRun);
end
% SAVE results
save(file_sim,'params','firings','PlotData');
fprintf('\nDone. Simulation saved in %s\n',file_sim);
toc
%%
%------------------------
% PLOT EVENTS
%------------------------
iTrial=1; % pick which trial to plot
dataload=load(file_sim); % load simulation results
data=dataload.PlotData{iTrial}; % membrane potentials
firings=dataload.firings{iTrial}; % spikes
Params=dataload.params; % parameters
Params.savedir=savedir;
plot_trial(data,firings,Params);
% figure 1 - rasterplot of all neurons in trial
% figure 2 - time course of membrane potential and PSC traces for E and I representative neurons
% figure 3 - time course of firing rate in clusters
% figure 4 - time course of stimuli (with CSgauss stimulus, the cue profile should be multiplied by a factor drawn from figure 5, one for each neuron
% figure 5 (with CSgauss stimulus only) - across-neurons distribution of cue peak values
