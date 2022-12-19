function [C]=C_assignment_MCMC(Z,A,W,C,S,rho_diag,rho_cc,par,alpha,type,method,gap_type,gap,opts)
% output: C (updated matrix of assignments)
    [~,J]=size(A{1});
    N=length(A);
    rho0p=par.rho0p;
    rho0n=par.rho0n;
    noc=size(Z,1);
    Ap=rho0p(1)*eye(noc)+rho0p(2)*(ones(noc)-eye(noc));
    An=rho0n(1)*eye(noc)+rho0n(2)*(ones(noc)-eye(noc)); 
    bs=round(J/50); % sample size at 2%
    accept=0;
    %for i=1:J % replace with s=1:bs for sampling
    for s=1:bs
        i=1+floor(J*rand);
        %disp(s);
        j=find(C(i,:)==1,1,'first'); %j: client with whom i was sitting
        l=find(Z(:,i)==1,1,'first'); %l: cluster to which i belonged
        g=1+floor(J*rand); % g: new client assignment, at random (coin toss)
        while g==j
            g=1+floor(J*rand);
        end
        B=C; % B: copy of C
        B(i,j)=0; % we move i from j to g
        B(i,g)=1;
        %B(g,i)=1;
        % We induce Z from B
        G=digraph(B);
        bins=conncomp(G, 'Type','weak');
        Z_copy=zeros(1,J);
        for i=1:J
            z=bins(1,i);
            Z_copy(z,i)=1;
        end
        k=find(Z_copy(:,i)==1,1,'first'); % now i belongs to k
        noc_copy=size(Z_copy,1);
        %Recompute some parameters for Z_copy (Ap, An, rho_diag, rho_cc)
        Ap_copy=rho0p(1)*eye(noc_copy)+rho0p(2)*(ones(noc_copy)-eye(noc_copy));
        An_copy=rho0n(1)*eye(noc_copy)+rho0n(2)*(ones(noc_copy)-eye(noc_copy)); 
        rho_diag_copy=zeros(noc_copy,N);
        if isfield(opts,'eta')
            for n=1:N
                rho_diag_copy(:,n)=diag(opts.eta(:,:,n));        
            end
        else
            for n=1:N
                switch method
                    case 'Binary'
                        rho_diag_copy(:,n)=betarnd(par.rho0p(1)*ones(noc_copy,1),par.rho0n(1)*ones(noc_copy,1));        
                    case 'Weighted'
                        rho_diag_copy(:,n)=gamrnd(par.rho0p(1)*ones(noc_copy,1),1/par.rho0n(1)*ones(noc_copy,1));        
                end
            end
        end
        if strcmp(gap_type,'same')
            rho_cc_copy=gap*minComDens(rho_diag_copy);
        else
            rho_cc_copy=minComDens(gap.*rho_diag_copy);      
        end
        % t-1: <i,j,l>; t: <i,g,k>; l->k
        [n_link,n_nonlink]=compute_links_between_clusters(A,Z,W,Ap,An,type,method); 
        [n_link_up,n_nonlink_up]=compute_links_between_clusters(A,Z_copy,W,Ap_copy,An_copy,type,method);
        [node_cl_p_bef,node_cl_np_bef]=compute_links_from_node_to_cluster(A,Z,W,type,method);
        [node_cl_p_aft,node_cl_np_aft]=compute_links_from_node_to_cluster(A,Z_copy,W,type,method);
        nilp=node_cl_p_bef(l,i,1);
        nilnp=node_cl_np_bef(l,i,1);
        nkip=node_cl_p_aft(k,i,1);
        nkinp=node_cl_np_aft(k,i,1);
        eta_ll=rho_diag(l);
        eta_kk=rho_diag_copy(k);
        v1=power(eta_ll,-nilp); 
        v2=power(1-eta_ll,-nilnp);
        v3=power(eta_kk,nkip);
        v4=power(1-eta_kk,nkinp);
        if j==i
            den=alpha;
        else
            den=S(i,j);
        end
        if g==i
            num=alpha;
        else
            num=S(i,g);
        end
        if den~=0
            v5=num/den;
        else
            v5=1;
        end
        v1=max(1e-323,v1);
        v2=max(1e-323,v2);
        v3=max(1e-323,v3);
        v4=max(1e-323,v4);
        v5=max(1e-323,v5);
        prod=v1*v2*v3*v4*v5;
        %
        noc=size(Z,1); % number of connected components
        for m=1:noc
            if (m~=l && m~=k)
                Nlmp=n_link(l,m,1);
                nimp=node_cl_p_bef(m,i,1);
                Nlmn=n_nonlink(l,m,1);
                nimn=node_cl_np_bef(m,i,1);
                xlm=rho_cc(l,m);
                %disp(xlm);
                fac1=max(0,Nlmp-nimp+rho0p(2));
                fac2=max(0,Nlmn-nimn+rho0n(2));
                num=betainc(xlm,fac1,fac2);
                fac1=max(0,Nlmp+rho0p(2));
                fac2=max(0,Nlmn+rho0n(2));
                den=betainc(xlm,fac1,fac2);
                if den~=0
                    quot=num/den;
                else
                    quot=1;
                end
                quot=max(1e-323,quot);
                %disp(quot);
                prod=prod*quot;
            end
        end
        %disp(prod);
        for m=1:noc_copy
            if (m~=k && m~=l)
                Nmkp=n_link_up(m,k);
                nmip=node_cl_p_aft(m,i,1);
                Nmkn=n_nonlink_up(m,k);
                nmin=node_cl_np_aft(m,i,1);
                xmk=rho_cc_copy(m,k);
                fac1=max(0,Nmkp+nmip+rho0p(2));
                fac2=max(0,Nmkn+nmin+rho0n(2));
                num=betainc(xmk,fac1,fac2);
                fac1=max(0,Nmkp+rho0p(2));
                fac2=max(0,Nmkn+rho0n(2));
                den=betainc(xmk,fac1,fac2);
                if den~=0
                    quot=num/den;
                else
                    quot=1;
                end
                quot=max(1e-323,quot);
                %disp(quot);
                prod=prod*quot;
            end
        end
        [Nlkp,nikp,nlip,Nlkn,nikn,nlin]=compute_lk(A{1},Z,Z_copy,k,l,i);
        %disp(['Nlkp: ' num2str(Nlkp) ' nikp: ' num2str(nikp) ' nlip: ' num2str(nlip) ' Nlkn: ' num2str(Nlkn) ' nikn: ' num2str(nikn) ' nlin: ' num2str(nlin)]);
        %xlk=rho_cc(l,k); pending!!!
        if noc_copy>noc
            xlk=rho_cc_copy(l,k);
        else
            xlk=rho_cc(l,k);
        end
        fac1=max(0,Nlkp-nikp+nlip+rho0p(2));
        fac2=max(0,Nlkn-nikn+nlin+rho0n(2));
        num=betainc(xlk,fac1,fac2);
        fac1=max(0,Nlkp+rho0p(2));
        fac2=max(0,Nlkn+rho0n(2));
        den=betainc(xlk,fac1,fac2);
        if den~=0
            quot=num/den;
        else
            quot=1;
        end
        tau=prod*quot;
        tau=min(1,tau);
        tau=max(0,tau);
        pp=cumsum(tau);
        if rand<pp
            %disp(' ');
            %disp(['Tau (accept): ' num2str(tau) ' | S(i,j): ' num2str(S(i,j)) ' | S(i,g): ' num2str(S(i,g)) ' | v5: ' num2str(v5)]);
            C=B;
            clear B;
            accept=accept+1;
        end
    end
    %disp([' Accepted ' num2str(accept) ' samples for C']);


function [n_link_up,n_nonlink_up] = compute_links_between_clusters(A,Z,W,Ap,An,type,method)
%links and non links from cluster i to cluster l conditioned on Z
    sumS=sum(Z,2);
    N=length(A);
    noc=length(sumS);
    SST=(sumS*sumS'-diag(sumS));
    n_link_up=zeros(noc,noc,N);
    n_nonlink_up=zeros(noc,noc,N);
    for n=1:N        
        if strcmp(type,'UnDirected')
            switch method
                case 'Binary'   
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';
                    n_link_up(:,:,n)=SASt+SASt';    
                    n_link_up(:,:,n)=n_link_up(:,:,n)-0.5*diag(diag(n_link_up(:,:,n)));
                    n_link_up(:,:,n)=n_link_up(:,:,n)+Ap;
                    n_nonlink_up(:,:,n)=SST-SASt-SWSt-SASt'-SWSt';
                    n_nonlink_up(:,:,n)=n_nonlink_up(:,:,n)-0.5*diag(diag(n_nonlink_up(:,:,n))); 
                    n_nonlink_up(:,:,n)=n_nonlink_up(:,:,n)+An;
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';
                    n_link_up(:,:,n)=SASt+SASt';    
                    n_link_up(:,:,n)=n_link_up(:,:,n)-0.5*diag(diag(n_link_up(:,:,n)));
                    n_nonlink_up(:,:,n)=SST-SWSt-SWSt';
                    n_nonlink_up(:,:,n)=n_nonlink_up(:,:,n)-0.5*diag(diag(n_nonlink_up(:,:,n)));                                            
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
            end  
            A{n}=A{n}+A{n}';
            W{n}=W{n}+W{n}';
        else                       
            switch method
                case 'Binary'                             
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link_up(:,:,n)=SASt+Ap;    
                    n_nonlink_up(:,:,n)=SST-SASt-SWSt+An;
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link_up(:,:,n)=SASt+Ap;    
                    n_nonlink_up(:,:,n)=SST-SWSt+An;                                        
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
            end                
        end        
    end
    
function [n_link_up,n_nonlink_up] = compute_links_from_node_to_cluster(A,Z,W,type,method)
%links and non links from node i to nodes in cluster l conditioned on Z
    [~,J]=size(A{1});
    N=length(A);
    sumS=sum(Z,2);
    noc=length(sumS);
    n_link_up=zeros(noc,J,N);
    n_nonlink_up=zeros(noc,J,N);
    sumS1=sum(Z,1);
    sumS2=sum(Z,2);
    SST=sumS2*sumS1;
    for n=1:N        
        if strcmp(type,'UnDirected')
            switch method
                case 'Binary'   
                    SASt=Z*A{n};
                    SWSt=Z*W{n};
                    n_link_up(:,:,n)=SASt;    
                    n_nonlink_up(:,:,n)=SST-SASt-SWSt; 
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n};
                    SWSt=Z*W{n};
                    n_link_up(:,:,n)=SASt;
                    n_nonlink_up(:,:,n)=SST-SWSt;
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
            end  
            A{n}=A{n}+A{n}';
            W{n}=W{n}+W{n}';
        else                        
            switch method
                case 'Binary'                             
                    SASt=Z*A{n};
                    SWSt=Z*W{n};                    
                    n_link_up(:,:,n)=SASt;    
                    n_nonlink_up(:,:,n)=SST-SASt-SWSt; 
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
                case 'Weighted'
                    SASt=Z*A{n};
                    SWSt=Z*W{n};                    
                    n_link_up(:,:,n)=SASt;    
                    n_nonlink_up(:,:,n)=SST-SWSt; 
                    for xx=1:noc
                        for xy=1:noc
                            if n_nonlink_up(xx,xy,n) < 0
                                n_nonlink_up(xx,xy,n) = 0;
                            end
                        end
                    end
            end                
        end        
    end

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

function [Nlkp,nikp,nlip,Nlkn,nikn,nlin]=compute_lk(A,Z,Z_copy,k,l,i)
nodes_in_l=find(Z(l,:)==1);
nodes_in_k=find(Z_copy(k,:)==1);
size_l=length(nodes_in_l);
size_k=length(nodes_in_k);
Total=size_l*size_k;
Nlkp=0;
for u=1:size_l
    n_l=nodes_in_l(u);
    for v=1:size_k
        n_k=nodes_in_k(v);
        %disp(['n_l: ' num2str(n_l) ' n_k: ' num2str(n_k)]);
        if A(n_l,n_k)==1
            Nlkp=Nlkp+1;
        end
    end
end
Nlkn=Total-Nlkp;
nikp=0;
nlip=0;
for u=1:size_l
    n_l=nodes_in_l(u);
    if A(n_l,i)==1
        nlip=nlip+1;
    end
end
for v=1:size_k
    n_k=nodes_in_k(v);
    if A(i,n_k)==1
        nikp=nikp+1;
    end
end
nikn=size_k-nikp;
nlin=size_l-nlip;

        

