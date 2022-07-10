% Ref [2] : NBS-based antenna selection in MIMO system

function  [capacityOfNBS, NBS_H_sel]=AS_NBS(Nr,Nt,Lr,SNR,H,fullAntenna)

% Record the largest channel capacity in all subsets when the number of selected antennas is Lr
capacityOfNBS=0;
selAntenna=fullAntenna;

% Correlation matrix
alpha=zeros(Nr,Nt);

% The degree of correlation between k and l
x=0;

%Antenna to be deleted
del=[];
for k=1:Nr
    for l=1:Nr
        if k<=l
            % The degree of correlation that does not need calculation or has been eliminated is recorded as -1
            alpha(k,l)=-1;
        else
            hk=H(k,:);
            hl=H(l,:);
            x=abs(dot(hk,hl));
            alpha(k,l)=x;
        end
    end
end
for m=Nr:-1:Lr+1
    [p, q]=find(max(max(alpha)));
    Xk=norm(H(p,:));
    Xl=norm(H(q,:));
    if Xk>=Xl
        alpha(q,:)=-1;
        alpha(:,q)=-1;
        del=[del,q];
    else
        alpha(p,:)=-1;
        alpha(:,p)=-1;
        del=[del,p];
    end
end
for n=1:length(del)
    x=del(n);
    
    % Delete antenna number x
    selAntenna=[selAntenna(1:x-1),selAntenna(x+1:end)];
end
NBS_H_sel=H(selAntenna,:);

% Channel capacity after selection
capacityOfNBS=log2(det(eye(Nt)+SNR/Nt*(NBS_H_sel'*NBS_H_sel)));
end

