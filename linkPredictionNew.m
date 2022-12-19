function [TP,TN,FP,FN,TPR,FPR]=linkPredictionNew(W,sample,A)
[i,j,~]=find(W);
iters=length(sample.iteration);
numLinks=length(i);
J=length(sample.Z{1});
b=1;
%numTrials=20;
TP=0;
FN=0;
for n=1:numLinks
    u=i(n);
    v=j(n);
    success=0;
    for m=1:iters
        Z_=sample.Z{m};
        eta_=sample.eta{m};
        gap_=sample.gap{m};
        u_cl=find(Z_(:,u));
        v_cl=find(Z_(:,v));
        if u_cl==v_cl
            %nu=eta_(u_cl,v_cl);
            %disp(nu);
            %if nu >=0.5
            success=success+1;
            %end
        else
            L_u=find(Z_(u_cl,:)==1);
            L_v=find(Z_(v_cl,:)==1);
            n_zi_zj_pos=0;
            n_zi_zj_neg=0;
            for x=1:length(L_u)
                for y=1:length(L_v)
                    x_i=L_u(x);
                    y_i=L_v(y);
                    if A(x_i,y_i)==1
                        n_zi_zj_pos=n_zi_zj_pos+1;
                    else
                        n_zi_zj_neg=n_zi_zj_neg+1;
                    end
                end
            end
            nu_1=eta_(u_cl,u_cl);
            nu_2=eta_(v_cl,v_cl);
            lim1=gap_*nu_1;
            lim2=gap_*nu_2;
            lim=min(lim1,lim2);
            alpha=n_zi_zj_pos+b+1;
            beta=n_zi_zj_neg+b;
            numerador=betainc(lim,alpha,beta);
            %disp('numerador');
            %disp(numerador);
            alpha=alpha-1;
            denominador=betainc(lim,alpha,beta);
            %disp('denominador');
            %disp(denominador);
            if numerador >= denominador
                success=success+1;
            end
        end
    end
    p=success/iters;
    %disp(p);
    if p >= 0.5
        TP=TP+1;
    end
end
FN=numLinks-TP;
TPR=(TP)/(TP+FN);
TN=0;
FP=0;
for n=1:numLinks
    neg_link=0;
    while neg_link==0
        u=ceil(rand*J);
        v=ceil(rand*J);
        if A(u,v)==0 && A(v,u)==0
            neg_link=1;
        end
    end
    success=0;
    for m=1:iters
        Z_=sample.Z{m};
        eta_=sample.eta{m};
        gap_=sample.gap{m};
        u_cl=find(Z_(:,u));
        v_cl=find(Z_(:,v));
        if u_cl==v_cl
            %nu=eta_(u_cl,v_cl);
            %if nu <= 0.5
            success=success+1;
            %end
        else
            L_u=find(Z_(u_cl,:)==1);
            L_v=find(Z_(v_cl,:)==1);
            n_zi_zj_pos=0;
            n_zi_zj_neg=0;
            for x=1:length(L_u)
                for y=1:length(L_v)
                    x_i=L_u(x);
                    y_i=L_v(y);
                    if A(x_i,y_i)==1
                        n_zi_zj_pos=n_zi_zj_pos+1;
                    else
                        n_zi_zj_neg=n_zi_zj_neg+1;
                    end
                end
            end
            nu_1=eta_(u_cl,u_cl);
            nu_2=eta_(v_cl,v_cl);
            lim1=gap_*nu_1;
            lim2=gap_*nu_2;
            lim=min(lim1,lim2);
            alpha=n_zi_zj_pos+b+1;
            beta=n_zi_zj_neg+b;
            numerador=betainc(lim,alpha,beta);
            %disp('numerador');
            %disp(numerador);
            alpha=alpha-1;
            denominador=betainc(lim,alpha,beta);
            %disp('denominador');
            %disp(denominador);
            if numerador <= denominador
                success=success+1;
            end
        end
    end
    p=success/iters;
    if p >= 0.5
        TN=TN+1;
    end
end
FP=numLinks-TN;
FPR=(FP)/(TN+FP);
