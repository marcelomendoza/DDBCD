function eta_norm=row_stochastic(eta)
eta_r=sum(eta);
R=length(eta);
eta_norm=zeros(R,R);
for i=1:R
    for j=1:R
        eta_norm(i,j)=eta(i,j)/eta_r(i);
    end
end