% This scripts generates a graph according to the Bayesian Community Detection (BCD) model and infer the
% parameters of the model using BCD and IRM

J=200; % Number of nodes, def: 250
nob=10; % Number of blobs

D2=generateDistanceMatrix(J,nob);

alpha=ceil(log(J)); % Parameter to C~DDCRP(alpha)
bp=[10 5];
%bn=[50 100]; %         eta_ll~Beta(bp(1),bn(1)), eta_lm~truncatedBeta(bp(2),bn(2))
bn=[5 50];
gap_prior=[5 20]; %      gap~Beta(gap_prior(1),gap_prior(1));
%type='UnDirected';
gap=0.1;

[A,Z_true,~,~,~,~,~,noc]=generateDDBCDGraph(J,D,nob,alpha,bp,bn,gap_prior,type,gap);
%disp(noc);


% Create Validation Data
pct_missing=2.5;
[W,class]=createValidationData(A,pct_missing,type);
%display(W);
% Run the DDBCD algorithm
opts.init_sample_iter=200; % Use 400 burn in iterations, discarded for MAP
opts.nsampleiter=100;      % Use 100 samples for MAP
%opts.type=type;
opts.alpha=ceil(log(J));
opts.gap=gap; % gamma
opts.dSstep=25; % Save every 50 samples for MAP

%DDBCD(A,D,W,noc,opts);
[L,cpu_time,NOC,eta_,gap_,Z,eta,gap,sample,West,predL,G]=DDBCD(A,D,W,noc,opts,bp,bn);
% Plot the results
[A_sorted,Z_sorted,eta_sorted,~]=sortGraphUnipartite(A,Z,eta);
plotSyntheticResults(A,West,Z_true,sample.MAP.Z,G,eta_sorted,L,NOC,eta_,gap_);

%% Run the IRM model
%[L_IRM,cpu_time_IRM,Z_IRM,eta_IRM,sample_IRM,West_IRM, par]=IRMUnipartite(A,W,noc,opts);
% Plot the results
%plotSyntheticResults(A,West,Z_true,sample_IRM.MAP.Z);

