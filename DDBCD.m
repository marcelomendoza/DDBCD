function [L,cpu_time,NOC,eta_,gap_,Z,eta,gap,sample,West,predL,G]=DDBCD(A,D,W,ANOC,opts,bp,bn)

% Non-parametric clustering of un-directed graphs based on MCMC
%
% Usage: 
%  [L,cpu_time,Z,eta,sample,West]=DDBCD(A,D,W,ANOC,opts)
%
% Input:
%   A       cell array of I x I  sparse adjacency matrix (can be triangular upper matrix)
%   D       distance matrix (square, JxJ)
%   W       cell array of I x I  sparse missing value indicator matrix (can be triangular upper matrix)
%   opts.
%           init_sample_iter initial number of sampling iterations  (i.e. burnin)
%           nsampleiter iterations to use for prediction
%           Z           I x noc matrix of community identities
%           dSstep      number of iterations between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           rho0p       1x2 vector of pseudo link counts wihtin and 
%                       between communities (default: [1 1])
%           rho0n       1x2 vector of pseudo non-link counts wihtin and 
%                       between communities and (default: [1 1])
%           gap_type    'same' or 'individual', (default: 'same')
%           constGap    make gap constant (default: false)
%           sameEta     use same value of eta for multiple graphs (default: false)

% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z           Estimated clustering assigment matrix
%   eta         Estimated relations
%   gap         Estimated Gap
%   sample      sampled parameters at each dSstep iterationCluster indicator matrix of size d x I where d<=noc before
%   West        Estimated link probabilities for missing links and non-links
%   predL       Predictive log-likelihood based on the use of the expected
%               value of eta rather than expected value of log(eta) and
%               log(1-eta) given as a sample average, i.e.
%               returns for Bernoulli: \sum_{ij\in W} A_ij*log(eta_ij)+(1-A_ij)*log(1-eta_ij)
%               returns for Poisson: \sum_{ij\in W} A_ij*log(eta_ij)-eta_ij
%

if nargin<5
    opts=struct;
end

% Initialize variables
if ~iscell(A)
    B=A;
    clear A;
    A{1}=B;
    clear B;
end
if ~iscell(W)
    B=W;
    clear W;
    W{1}=B;
    clear B;
end
par.N=length(A);
[I,J]=size(A{1}); % If A is an adjacency matrix, I=J
type=mgetopt(opts,'type','UnDirected');
method=mgetopt(opts,'method','Binary');
gap_type=mgetopt(opts,'gap_type','same');
par.constGap=mgetopt(opts,'constGap',false);
par.sameEta=mgetopt(opts,'sameEta',false);
alpha=opts.alpha;
Iw=cell(1,par.N);
Jw=cell(1,par.N);
West=cell(1,par.N);
sumA=0;
nnzA=0;
switch method
    case {'Binary'}
        par.rho0p=mgetopt(opts,'rho0p',[bp(1) bp(2)]); % pseudo counts of links between clusters
        par.rho0n=mgetopt(opts,'rho0n',[bn(1) bn(2)]); % pseudo counts of non-links between clusters        
    case {'Weighted'}
        par.rho0p=mgetopt(opts,'rho0p',[1 1]); % pseudo counts of links between clusters
        par.rho0n=mgetopt(opts,'rho0n',[1 1]); % pseudo counts of non-links between clusters                        
end


for n=1:par.N    
    if strcmp(type,'UnDirected')
        A{n}=triu(A{n});
        W{n}=triu(W{n});
    end    
    W{n}=logical(W{n});
    [~,~,vv1{n}]=find(W{n}.*A{n}+W{n});
    vv1{n}=vv1{n}-1;
    switch method
        case {'Binary'}
            A{n}=logical(A{n}-A{n}.*W{n}); % Remove missing links
        case {'Weighted'}
            A{n}=A{n}-A{n}.*W{n}; % Remove missing links
    end
    [Iw{n},Jw{n}]=find(W{n});
    West{n}=sparse(I,J);    
    sumA=sumA+sum(sum(A{n}));
    nnzA=nnzA+nnz(A{n});
    N=par.N;
    if par.sameEta
        if n==1
           At=A{n};
           Wt=W{n};
        else
            At=At+A{n};
            Wt=Wt+W{n};
        end
        if n==par.N
           clear A W
           A{1}=At;
           W{1}=Wt;
        end
        N=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate P (exponential decay)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1./(ANOC+1);
S=exp(-D./a);
s=size(S);
idx=1:s(1)+1:s(1)*s(2);
S(idx)=S(idx)*alpha;
P=zeros(J,J);
for i=1:J
    row=S(i,:)./sum(S(i,:));
    P(i,:)=row;
end
%disp('P ready');

%%%%%%%%%%%%%%%%%%%%%%%
% Generate assignments%
%%%%%%%%%%%%%%%%%%%%%%%
% We allocate data points conditioned on distance measurements

C=zeros(J,J);
for i=1:J % arrivals in order
   p=P(i,:);
   pp=cumsum(p);
   ind=find(rand<pp,1,'first'); % assignment according to pp
   C(i,ind)=1;  
   %C(ind,i)=1;
end
%disp('C ready');

%%%%%%%%%%%%%%%%%%%%
%Generate clusters %
%%%%%%%%%%%%%%%%%%%%
% We induce Z from C
G=digraph(C);
bins=conncomp(G, 'Type','weak');
Z=zeros(1,J);
for i=1:J
    z=bins(1,i);
    Z(z,i)=1;
end
noc=size(Z,1); % number of connected components
%disp('Z ready');

% Set remaining parameters

if strcmp(gap_type,'same')
    gap=mgetopt(opts,'gap',1);
else
    gap=mgetopt(opts,'gap',ones(noc,N));
end

predL=zeros(1,par.N);
rho_diag=zeros(noc,N);
if isfield(opts,'eta')
    for n=1:N
        rho_diag(:,n)=diag(opts.eta(:,:,n));        
    end
else
    for n=1:N
        switch method
            case 'Binary'
                rho_diag(:,n)=betarnd(par.rho0p(1)*ones(noc,1),par.rho0n(1)*ones(noc,1));        
            case 'Weighted'
                rho_diag(:,n)=gamrnd(par.rho0p(1)*ones(noc,1),1/par.rho0n(1)*ones(noc,1));        
        end
        
    end
end
%display(rho_diag);

alpha=mgetopt(opts,'alpha',log(J));
verbose=mgetopt(opts,'verbose',1);
dSstep=mgetopt(opts,'dSstep',25); % pseudo counts of links between clusters
init_sample_iter=mgetopt(opts,'init_sample_iter',25);
nsampleiter=mgetopt(opts,'nsampleiter',50);
maxiter=nsampleiter+init_sample_iter;

sample=struct;
L=nan(1,maxiter);
cpu_time=L;
NOC=L;
eta_=L;
gap_=L;
sstep=0;
westiter=0;
    
iter=0;
if verbose % Display algorithm    
    disp(['Non-parametric clustering based on the Bayesian Community Detection using ' method ' method for graphs']);
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','logP','dlogP/|logP|','noc','time','\n');
    dline = sprintf('+');
    disp(dline);
    fprintf(dheader);
    disp(dline);
end

Q=-inf;
Qbest=Q;

flag=1;
while iter<maxiter
    iter=iter+1;
    if (noc==ANOC && flag==1)
        maxiter=iter+nsampleiter;
        flag=0;
        disp('--------ANOC detected--------');
    end
    if (iter==49500) % max number of iterations
        maxiter=iter+nsampleiter;
        flag=0;
        disp('--------ANOC detected--------');
    end
    tic;
    Qold=Q;               
   

    if strcmp(gap_type,'same')
    rho_cc=gap*minComDens(rho_diag);
    else
    rho_cc=minComDens(gap.*rho_diag);      
    end

    % Gibbs sampling of Z in random order        
    %[Z,rho_diag,n_link,n_nonlink,gap]=Gibbs_sample_SBGC(gap,Z,A,W,rho_diag,par,alpha,randperm(J),type,method,gap_type);
    C=C_assignment_MCMC(Z,A,W,C,S,rho_diag,rho_cc,par,alpha,type,method,gap_type,gap,opts);     
    %disp('C_MCMC assignment ready');

    % We induce Z from C; Z(C)
    G=digraph(C);
    bins=conncomp(G, 'Type','weak');
    Z=zeros(1,J);
    for i=1:J
        z=bins(1,i);
        Z(z,i)=1;
    end
    noc=size(Z,1);
    sumS=sum(Z,2);
    N=length(A); % number of connected components
    %disp('Z inference ready');
    
    % We compute rho_diag, rho_cc (gap), n_link and n_nonlink from Z
    rho0p=par.rho0p;
    rho0n=par.rho0n;
    Ap=rho0p(1)*eye(noc)+rho0p(2)*(ones(noc)-eye(noc));
    An=rho0n(1)*eye(noc)+rho0n(2)*(ones(noc)-eye(noc));
    rho_diag=zeros(noc,N);
    if isfield(opts,'eta')
        for n=1:N
            rho_diag(:,n)=diag(opts.eta(:,:,n));        
        end
    else
        for n=1:N
            switch method
                case 'Binary'
                    rho_diag(:,n)=betarnd(par.rho0p(1)*ones(noc,1),par.rho0n(1)*ones(noc,1));        
                case 'Weighted'
                    rho_diag(:,n)=gamrnd(par.rho0p(1)*ones(noc,1),1/par.rho0n(1)*ones(noc,1));        
            end
        end
    end
    
    SST=(sumS*sumS'-diag(sumS));
    %display(sumS);
    n_link=zeros(noc,noc,N);
    n_nonlink=zeros(noc,noc,N);
    for n=1:N        
        if strcmp(type,'UnDirected')
            switch method
                case 'Binary' 
                    %display(Z);
                    %display(A{1});
                    %display(W{1});
                    SASt=Z*A{n}*Z';
                    %display(SASt);
                    SWSt=Z*W{n}*Z';
                    %display(SWSt);
                    n_link(:,:,n)=SASt+SASt';    
                    %display(n_link);
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    %display(n_link);
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    %display(n_link);
                    %display(SST);
                    n_nonlink(:,:,n)=SST-SASt-SWSt-SASt'-SWSt';
                    %display(n_nonlink);
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    %display(n_nonlink);
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;
                    %display(n_nonlink);
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink(xx,xy,n) < 0
                                n_nonlink(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';
                    n_link(:,:,n)=SASt+SASt';    
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    n_nonlink(:,:,n)=SST-SWSt-SWSt';
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink(xx,xy,n) < 0
                                n_nonlink(xx,xy,n) = 0;
                            end
                        end
                    end
            end 
            A{n}=logical(A{n}+A{n}');
            %display(A{1});
        else                       
            switch method
                case 'Binary'                             
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link(:,:,n)=SASt+Ap;    
                    n_nonlink(:,:,n)=SST-SASt-SWSt+An; 
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink(xx,xy,n) < 0
                                n_nonlink(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link(:,:,n)=SASt+Ap;    
                    n_nonlink(:,:,n)=SST-SWSt+An;  
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink(xx,xy,n) < 0
                                n_nonlink(xx,xy,n) = 0;
                            end
                        end
                    end
            end                
        end        
    end
    %display(Z);
    %display(n_link);
    %display(n_nonlink);
    % MH sample rho_diag    
    rho_diag=MHsampleRho_diag(gap,n_link,n_nonlink,rho_diag,type,method,par);
    
    % MH sample gap
    [logP_A,logP_S]=evalLikelihood(gap,n_link,n_nonlink,rho_diag,par,alpha,Z,type,method);        
    if ~par.constGap
        [gap,logP_A]=MHsampleGap(gap,gap_type,n_link,n_nonlink,rho_diag,type,method,par.rho0p,par.rho0n,logP_A);        
    end
        
    % form expected value of rho_noise by numeric integration
    rho_noise=zeros(noc,noc,N);
    if strcmp(gap_type,'same')
        rho_cc=gap*minComDens(rho_diag);      
    else
        rho_cc=minComDens(gap.*rho_diag);      
    end
    if strcmp(type,'UnDirected')
        ii=find(triu(ones(noc),1));
    else
        ii=find(ones(noc)-eye(noc));
    end
          
    for n=1:N
        n_link_t=n_link(:,:,n);
        n_nonlink_t=n_nonlink(:,:,n);
        rho_cc_t=rho_cc(:,:,n);
        rho_noise_t=zeros(noc);
        switch method            
            case 'Binary'
                %rho_noise=(beta(n_link+1,n_nonlink).*betainc(rho_cc,n_link+1,n_nonlink))./(beta(n_link,n_nonlink).*betainc(rho_cc,n_link,n_nonlink));                
                rho_noise_t(ii)=n_link_t(ii)./(n_link_t(ii)+n_nonlink_t(ii)).*betainc(rho_cc_t(ii),n_link_t(ii)+1,n_nonlink_t(ii))./(betainc(rho_cc_t(ii),n_link_t(ii),n_nonlink_t(ii))+1e-320);            
                idx=find(rho_noise_t(ii)==0);
                % if density too low use numerical integration
                rho_noise_t(ii(idx))=my_noise_est_beta(rho_cc_t(ii(idx)),n_link_t(ii(idx)),n_nonlink_t(ii(idx)));                                     
            case 'Weighted'
                %rho_noise=(gamma(n_link+1)./(n_nonlink.^(n_link+1)).*gammainc(rho_cc.*n_nonlink,n_link+1))./(gamma(n_link)./(n_nonlink.^n_link).*gammainc(rho_cc.*n_nonlink,n_link));            
                rho_noise_t(ii)=(n_link_t(ii)./n_nonlink_t(ii)).*gammainc(rho_cc_t(ii).*n_nonlink_t(ii),n_link_t(ii)+1)./(gammainc(rho_cc_t(ii).*n_nonlink_t(ii),n_link_t(ii))+1e-320);                            
                idx=find(rho_noise_t(ii)==0);
                % if density is too low use numerical integration
                rho_noise_t(ii(idx))=my_noise_est_gamma(rho_cc_t(ii(idx)),n_link_t(ii(idx)),n_nonlink_t(ii(idx)));                                 
        end
        rho_noise(:,:,n)=rho_noise_t;
    end
           
    eta=zeros(noc,noc,N);
    for n=1:length(A) 
         if strcmp(type,'UnDirected')
             eta(:,:,n)=rho_noise(:,:,n)+rho_noise(:,:,n)'+diag(rho_diag(:,n));         
         else
             eta(:,:,n)=rho_noise(:,:,n)+diag(rho_diag(:,n));         
         end
    end              
    
    Q=logP_A+logP_S;
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;  
    NOC(iter)=noc;
    eta_(iter)=eta(1,1);
    gap_(iter)=gap(1,1);
    if mod(iter,dSstep)==0 && flag==0
        sstep=sstep+1;
        sample.iteration(sstep)=iter;
        sample.Z{sstep}=Z;        
        sample.eta{sstep}=eta;  
        sample.gap{sstep}=gap;  
    end
    if Q>Qbest
        Qbest=Q;
        sample.MAP.L=Q;
        sample.MAP.iteration=iter;
        sample.MAP.Z=Z;
        sample.MAP.eta=eta;
        sample.MAP.gap=gap;
    end
    if rem(iter,1)==0 && verbose && mod(iter,dSstep)==0      
        fprintf('It: %4.0f | Q: %8.4e | dQ/abs(Q): %8.4e | noc: %4.0f | t_iter: %8.4f \n',iter,Q,dQ/abs(Q),noc,t_iter);
    end
        
    % Estimate missing link probability of sample    
    if  flag==0 
        westiter=westiter+1;        
        if iter==maxiter-nsampleiter+1
            %disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iterations']);   
        end
        step=10000;        
        for n=1:par.N            
            val=zeros(1,length(Iw{n}));
            for k=1:ceil((length(Iw{n})/step))
                ind=(k-1)*step+1:min([k*step, length(Iw{n})]);
                if par.sameEta
                    val(ind)=sum(Z(:,Iw{n}(ind)).*(eta*Z(:,Jw{n}(ind))))+eps;
                else
                    val(ind)=sum(Z(:,Iw{n}(ind)).*(eta(:,:,n)*Z(:,Jw{n}(ind))))+eps;
                end
            end
            West{n}=West{n}+sparse(Iw{n},Jw{n},val,I,J)/nsampleiter;                
            if strcmp(method,'Weighted')
                predL(n)=predL(n)+sum(vv1{n}.*log(val')-val')/(nsampleiter*length(val));
            else                
                predL(n)=predL(n)+sum(vv1{n}.*log(val')+(1-vv1{n}).*log(1-val'))/(nsampleiter*length(val));
            end
        end
    end
     %   figure(1); subplot(1,2,1); plot(L); subplot(1,2,2); imagesc(rho_diag); colorbar;
end

% sort the communities
[~,ind]=sort(sum(Z,2),'descend');
Z=Z(ind,:);
eta=eta(ind,:,:);
eta=eta(:,ind,:);
if ~strcmp(gap_type,'same')
   gap=gap(ind); 
end
if verbose   
  fprintf('%12.0f | %12.4e | %12.4e | %12.0f |%12.4f \n',iter,Q,dQ/abs(Q),noc,t_iter);
  disp('');
end
%}
            
% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset'
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise
            error('Wrong option: %s.', cmd);
    end
end

% -------------------------------------------------------------------------
 function rho_diag=MHsampleRho_diag(gap,n_link,n_nonlink,rho_diag,type,method,par)
        
     rho0p=par.rho0p;
     rho0n=par.rho0n;
     if strcmp(method,'Binary')
        clusterFun=@my_betaub;
        logLike=@lnbetalike;
        drawsample=@betarnd;
     else
        clusterFun=@my_gammaub;
        logLike=@lngammalike;
        drawsample=@gamrnd;
     end
    M=100;
    N=size(rho_diag,2); 
    noc=size(n_link,1);
    accept=0;
    sigma_sq=0.01;
    for rep=1:10  % sample for diagonal elements of eta many times since this is inexpensive
        for n=1:N
            rho_diag_t=rho_diag(:,n);            
            n_link_t=n_link(:,:,n);
            n_nonlink_t=n_nonlink(:,:,n);
            if size(gap,2)>1
               gap_t=gap(:,n); 
            else
               gap_t=gap;
            end
            for d=randperm(noc)
                if strcmp(type,'UnDirected')
                   B=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                else
                   B1=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                   B2=clusterFun(gap_t,n_link_t(d,:)',n_nonlink_t(d,:)',rho_diag_t,rho0p,rho0n,d);
                   B=B1+B2;
                   B(d)=B1(d);
                end
                if rem(rep,3)==1
                    rho_diag_t(d)=drawsample(1,1);
                    deltalogQ=0;
                elseif rem(rep,3)==2
                    a=n_link(d,d);
                    b=n_nonlink(d,d);
                    if strcmp(method,'Weighted')                        
                        rho_diag_t(d)=drawsample(a,1/b);
                    else
                        rho_diag_t(d)=drawsample(a,b);
                    end                    
                    deltalogQ=logLike(rho_diag(d,n),a,b)-logLike(rho_diag_t(d),a,b);                        
                elseif rem(rep,3)==0                                        
                    if strcmp(method,'Weighted')                        
                        mu=rho_diag(d,n);
                        a=mu^2/sigma_sq;
                        b=sigma_sq/mu;            
                        rho_diag_t(d)=drawsample(a,b);
                        mu_new=rho_diag_t(d);
                        a_new=mu_new^2/sigma_sq;
                        b_new=sigma_sq/mu_new;            
                    else                        
                        a=M*rho_diag(d,n);
                        b=M-a;
                        rho_diag_t(d)=drawsample(a,b);
                        a_new=M*rho_diag_t(d);
                        b_new=M-a_new;
                    end                    
                    deltalogQ=logLike(rho_diag(d,n),a_new,b_new)-logLike(rho_diag_t(d),a,b);                        
                end          
                
                if strcmp(type,'UnDirected')
                   Bnew=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                else
                   B1=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                   B2=clusterFun(gap_t,n_link_t(d,:)',n_nonlink_t(d,:)',rho_diag_t,rho0p,rho0n,d);
                   Bnew=B1+B2;
                   Bnew(d)=B1(d); % subtract contribution from d to d which has been included twice (fixed 12 August 2013)
                end
                if rand<exp(sum(Bnew)-sum(B)+deltalogQ) % notice log proposal = 0 hence ignored in MH move
                    accept=accept+1;
                    rho_diag(d,n)=rho_diag_t(d);
                else
                    rho_diag_t(d)=rho_diag(d,n);
                end                
            end
        end
    end
    %disp([' Accepted ' num2str(accept) ' samples for eta']);
        

% -------------------------------------------------------------------------  
function [logP_A,logP_S]=evalLikelihood(gap,n_link,n_nonlink,rho_diag,par,alpha,Z,type,method)
    
    rho0p=par.rho0p;
    rho0n=par.rho0n;
    N=size(n_link,3);
    [noc, J]=size(Z);    
    logP_A=0;            
    if strcmp(method,'Binary')
        clusterFun=@my_betaub;        
    else
        clusterFun=@my_gammaub;        
    end
    for n=1:N                
        if strcmp(type,'UnDirected')            
            if size(gap,2)>1
                logP_A=logP_A+sum(sum(triu(clusterFun(gap(:,n),n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n))));
            else
                logP_A=logP_A+sum(sum(triu(clusterFun(gap,n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n))));
            end
        else
            if size(gap,2)>0
                 logP_A=logP_A+sum(sum(clusterFun(gap(:,n),n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n)));
            else
                logP_A=logP_A+sum(sum(clusterFun(gap,n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n)));
            end
        end        
    end
    sumS=sum(Z,2);        
    logP_S=noc*log(alpha)+sum(gammaln(full(sumS)))-gammaln(J+alpha)+gammaln(alpha);    
    
 % -------------------------------------------------------------------------
 function [gap,logP]=MHsampleGap(gap,gap_type,n_link,n_nonlink,rho_diag,type,method,rho0p,rho0n,logP)
        
    if strcmp(method,'Binary')
        clusterFun=@my_betaub;        
    else
        clusterFun=@my_gammaub;        
    end    
    accept=0;    
    [noc,~, N]=size(n_link);
    for rep=1:10  % sample for gap many times since this is inexpensive                                
        M=ceil(250*rand);
        if strcmp(gap_type,'same')
            if rem(rep,2)==1
                gap_t=betarnd(1,1);
                deltalogQ=0;
            elseif rep==2   
                a=M*gap;
                b=M-a;                                    
                gap_t=betarnd(a,b);                           
                a_new=M*gap_t;
                b_new=M-a_new;                                    
                deltalogQ=lnbetalike(gap_t,a,b)-lnbetalike(gap,a_new,b_new);                                        
            end                                
           Bnew=clusterFun(gap_t,n_link,n_nonlink,rho_diag,rho0p,rho0n);
           logPnew=0;
           for n=1:N
               if strcmp(type,'UnDirected')
                   logPnew=logPnew+sum(sum(triu(Bnew(:,:,n))));
               else
                   logPnew=logPnew+sum(sum(Bnew(:,:,n)));
               end
           end
            if rand<exp(logPnew-logP+deltalogQ) % notice log proposal = 0 hence ignored in MH move
               accept=accept+1;
               gap=gap_t;
               logP=logPnew;
            end                           
        else
            for n=1:N
                gap_t=gap(:,n);
                for d=randperm(noc)            
                   if strcmp(type,'UnDirected')
                       B=clusterFun(gap,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                    else
                       B1=clusterFun(gap(:,n),n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                       B2=clusterFun(gap(:,n),n_link(d,:,n)',n_nonlink(d,:,n)',rho_diag(:,n),rho0p,rho0n,d);
                       B=B1+B2;
                       B(d)=B1(d);
                    end
                    if rem(rep,2)==1
                        gap_t(d)=betarnd(1,1);
                        deltalogQ=0;
                    else                        
                        a=M*gap(d,n);
                        b=M-a;                                    
                        gap_t(d)=betarnd(a,b);                           
                        a_new=M*gap_t(d);
                        b_new=M-a_new;                                    
                        deltalogQ=lnbetalike(gap_t(d),a,b)-lnbetalike(gap_t(d),a_new,b_new);                        
                    end          

                    if strcmp(type,'UnDirected')
                       Bnew=clusterFun(gap_t,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                    else
                       B1=clusterFun(gap_t,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                       B2=clusterFun(gap_t,n_link(d,:,n)',n_nonlink(d,:,n)',rho_diag(:,n),rho0p,rho0n,d);
                       Bnew=B1+B2;
                       Bnew(d)=B1(d);
                    end
                    if rand<exp(sum(Bnew)-sum(B)+deltalogQ) % notice log proposal = 0 hence ignored in MH move
                        accept=accept+1;
                        gap(d,n)=gap_t(d);
                    else
                        gap_t(d)=gap(d,n);
                    end                
                end 
            end
        end
    end
    gap(gap<eps)=eps; % Ensure numeric stability
    %disp([' Accepted ' num2str(accept) ' samples for gap']);
   
   

% --------------------------------------
 function B=my_betaub(gap,n_link,n_nonlink,rho_diag,rho0p,rho0n,dd)
        % n_link        D x DD x N matrix
        % n_nonlink     D x DD x N matrix
        % rho_diag      D x N matrix 
        % dd             indices corresponding to DD
                
        [D,DD,N]=size(n_link);                             
        beta_diag_const=betaln(rho0p(1),rho0n(1));
        
        if nargin<7
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag);           
                E=ones(size(rho_diag));            
                BetaIncp=betainc(gap*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            else
                rho_min=minComDens(gap.*rho_diag);           
                E=ones(size(rho_diag));           
                BetaIncp=betainc(gap.*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            end
            
            prior_min=minComDens(prior_diag);
            BetaInc=betainc(rho_min(1:size(n_link,1),:,:),n_link,n_nonlink);
            BetaInc(BetaInc==0)=1e-323;
            logBetaInc=log(BetaInc);
            B = betaln(n_link,n_nonlink)+logBetaInc-prior_min(1:size(n_link,1),:,:);        
            for n=1:N
                if D==DD
                    dn_link=diag(n_link(:,:,n));
                    dn_nonlink=diag(n_nonlink(:,:,n));
                    ii=find(diag(B(:,:,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(:,:,n)=B(:,:,n)-diag(diag(B(:,:,n)))+diag(log(rho_diag(:,n)).*(dn_link-1)+log(1-rho_diag(:,n)).*(dn_nonlink-1)-beta_diag_const);           
                else
                    dn_link=diag(n_link(1:D,1:D,n));
                    dn_nonlink=diag(n_nonlink(1:D,1:D,n));
                    ii=find(diag(B(1:D,1:D,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(1:D,1:D,n)=B(1:D,1:D,n)-diag(diag(B(1:D,1:D,n)))+diag(log(rho_diag(1:D,n)).*(dn_link-1)+log(1-rho_diag(1:D,n)).*(dn_nonlink-1)-beta_diag_const);           
                end
            end                   
        else     
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag,dd);   
                E=ones(size(rho_diag));
                BetaIncp=betainc(gap*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            else
                rho_min=minComDens(gap.*rho_diag,dd);          
                E=ones(size(rho_diag));
                BetaIncp=betainc(gap.*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;                
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            end
            
            prior_min=minComDens(prior_diag,dd);   
            BetaInc=betainc(rho_min(1:size(n_link,1),:,:),n_link,n_nonlink);
            BetaInc(BetaInc==0)=1e-323;
            logBetaInc=log(BetaInc);
            B = betaln(n_link,n_nonlink)+logBetaInc-prior_min(1:size(n_link,1),:,:);                    
            if dd<=D                            
                for n=1:N
                    dn_link=diag(n_link(dd,:,n));
                    dn_nonlink=diag(n_nonlink(dd,:,n));
                    ii=find(diag(B(dd,:,n))==-Inf);
                    for t=1:length(ii)
                        B(dd(ii(t)),ii(t),n)=0;
                    end
                    B(dd,:,n)=B(dd,:,n)-diag(diag(B(dd,:,n)))+diag(log(rho_diag(dd,n)).*(dn_link-1)+log(1-rho_diag(dd,n)).*(dn_nonlink-1)-beta_diag_const);           
                end
            end            
        end
        
% --------------------------------------
 function B=my_gammaub(gap,n_link,n_nonlink,rho_diag,rho0p,rho0n,dd)
        % n_link        D x DD x N matrix
        % n_nonlink     D x DD x N matrix
        % rho_diag      D x N matrix 
        % dd             indices corresponding to DD
        
        [D,DD,N]=size(n_link);             
        gamma_diag_const=rho0p(1).*log(rho0n(1))-gammaln(rho0p(1));
        if nargin<7
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            else
                rho_min=minComDens(gap.*rho_diag);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap.*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            end            
            prior_min=minComDens(prior_diag);            
            GamInc=gammainc(rho_min(1:size(n_link,1),:,:).*n_nonlink,n_link);
            GamInc(GamInc==0)=1e-323;
            logGamInc=log(GamInc);            
            B = -n_link.*log(n_nonlink)+gammaln(n_link)+logGamInc-prior_min(1:size(n_link,1),:,:);        
            for n=1:N
                if D==DD
                    dn_link=diag(n_link(:,:,n));
                    dn_nonlink=diag(n_nonlink(:,:,n));
                    ii=find(diag(B(:,:,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t),n)=0;
                    end
                    B(:,:,n)=B(:,:,n)-diag(diag(B(:,:,n)))+diag(log(rho_diag(:,n)).*(dn_link-1)-rho_diag(:,n).*dn_nonlink-gamma_diag_const);           
                else
                    dn_link=diag(n_link(1:D,1:D,n));
                    dn_nonlink=diag(n_nonlink(1:D,1:D,n));
                    ii=find(diag(B(1:D,1:D,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(1:D,1:D,n)=B(1:D,1:D,n)-diag(diag(B(1:D,1:D,n)))+diag(log(rho_diag(1:D,n)).*(dn_link-1)-rho_diag(1:D,n).*dn_nonlink-gamma_diag_const);           
                end
            end                   
        else
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag,dd);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            else
                rho_min=minComDens(gap.*rho_diag,dd);                    
                E=ones(size(rho_diag));                
                GamIncp=gammainc(gap.*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            end            
            
            prior_min=minComDens(prior_diag,dd);     
            GamInc=gammainc(rho_min(1:size(n_link,1),:,:).*n_nonlink,n_link);
            GamInc(GamInc==0)=1e-323;
            logGamInc=log(GamInc);
            B = -n_link.*log(n_nonlink)+gammaln(n_link)+logGamInc-prior_min(1:size(n_link,1),:,:);
            if dd<=D  
                for n=1:N
                    dn_link=diag(n_link(dd,:,n));
                    dn_nonlink=diag(n_nonlink(dd,:,n));
                    ii=find(diag(B(dd,:,n))==-Inf);
                    for t=1:length(ii)
                        B(dd(ii(t)),ii(t),n)=0;
                    end
                    B(dd,:,n)=B(dd,:,n)-diag(diag(B(dd,:,n)))+diag(log(rho_diag(dd,n)).*(dn_link-1)-rho_diag(dd,n).*dn_nonlink-gamma_diag_const);           
                end
            end
        end

% ----------------------------------------------------
function rho_min=minComDens(rho_diag,d)
[noc,N]=size(rho_diag);
if nargin<2
    ed=ones(1,noc);
    rho_cc=zeros(noc,noc,N,2);
else
    ed=ones(1,length(d));
    enoc=ones(noc,1);
    rho_cc=zeros(noc,length(d),N,2);
end
rho_c=permute(rho_diag(:,:,ed),[1 3 2]);
rho_cc(:,:,:,1)=rho_c;
if nargin<2
    rho_cc(:,:,:,2)=permute(rho_c,[2 1 3]);   
else
    rho_cc(:,:,:,2)=permute(rho_diag(d,:,enoc),[3 1 2]);
end
rho_min=min(rho_cc,[],4);

% ----------------------------------------------------
function logL=lnbetalike(x,a,b)    
    logL=gammaln(a+b)-gammaln(a)-gammaln(b)+(a-1).*log(x+1e-323)+(b-1).*log(1-x+1e-323);

% ----------------------------------------------------
function logL=lngammalike(x,a,b)    
    logL=a.*log(b)-gammaln(a)+(a-1).*log(x)-b.*x;
            
%-----------------------------------------------------
    function noise_est=my_noise_est_beta(gap,n_link,n_nonlink)
        noise_est=nan(size(n_link));
        steps=10000;
        for t=1:length(n_link)
            x=linspace(1e-320,gap(t),steps);
            q=(n_nonlink(t)-1)*log(1-x);
            logP1=(n_link(t)+1-1)*log(x)+q;
            logP2=(n_link(t)-1).*log(x)+q;        
            maxlogP1=max(logP1);
            logP1=logP1-maxlogP1;
            logP2=logP2-maxlogP1;
            noise_est(t)=trapz(exp(logP1))/trapz(exp(logP2));
        end
        
%-----------------------------------------------------
    function noise_est=my_noise_est_gamma(gap,n_link,n_nonlink)
        noise_est=nan(size(n_link));
        steps=10000;
        for t=1:length(n_link)
            x=linspace(1e-320,gap(t),steps);        
            logP1=(n_link(t)+1-1)*log(x)-n_nonlink(t)*x;
            logP2=(n_link(t)-1).*log(x)-n_nonlink(t)*x;        
            maxlogP1=max(logP1);
            logP1=logP1-maxlogP1;
            logP2=logP2-maxlogP1;
            noise_est(t)=trapz(exp(logP1))/trapz(exp(logP2));
        end
               