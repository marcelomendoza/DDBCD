% {cornell, texas, washington, wisconsin, citeseer, cora}
eta=2; %2 for small networks, 3 or 4 for large networks
[IDs,IDMap,D,dim]=readContentFilewebkb('Datasets/webkb/cora.content', eta);
A=readGraphFilewebkb('Datasets/webkb/cora.cites',IDMap);
disp('Data processed');
J=length(A);
%alpha=2;
%alpha=ceil(log(J)); % Parameter to C~DDCRP(alpha)%
%bp=[10 5];
bp=[1 1];
%bn=[5 50];
bn=[1 1];
%gap_prior=[5 20]; %      gap~Beta(gap_prior(1),gap_prior(1));
%gap_prior=[1 1];
type='UnDirected';
gap=0.01;%0.01?
noc=7;

% Create Validation Data
pct_missing=10;
[W,class]=createValidationData(A,pct_missing,type);
disp('Validation data ready')

% Run the DDBCD algorithm
opts.init_sample_iter=50000; % Use 150 burn in iterations, discarded for MAP (200 for citeseer)
opts.nsampleiter=500;      % Use 100 samples for MAP
opts.type=type;
opts.alpha=ceil(log(J));
%opts.alpha=2;
opts.gap=gap; % gamma
opts.dSstep=5; % Save every 5 samples for MAP

%DDBCD(A,D,W,noc,opts);
[L,cpu_time,NOC,eta_,gap_,Z,eta,gap,sample,West,predL,G]=DDBCD(A,D,W,noc,opts,bp,bn);

% Plot the results
[A_sorted,Z_sorted,eta_sorted,~]=sortGraphUnipartite(A,Z,eta);
%plotSyntheticResults_extended(A,West,sample.MAP.Z,G,eta_sorted,L,NOC,eta_,gap_);

% DDBCD Evaluation
[TP,TN,FP,FN,TPR,FPR]=linkPredictionNew(W,sample,A);
disp('-------Results-------');
disp(TP);
disp(TN);
disp(FP);
disp(FN);
disp(TPR);
disp(FPR);
AUC=0.5-(FPR/2)+(TPR/2);
disp(AUC);
ACC=(TP+TN)/(TP+TN+FP+FN);
disp(ACC);




