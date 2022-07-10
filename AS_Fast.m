% Ref [3] : Fast antenna subset selection in MIMO System

function  [capacityOfSelected, Fast_H_sel]=AS_Fast(Nr,Nt,Lr,SNR,H,fullAntenna)

if(Lr==Nr)
    capacityOfSelected=log2(det(eye(Nt)+SNR/Nt*(H'*H)));
else
    
    B=eye(Nt,Nt);
    Alpha=[]; % record every alpha (empty array)
    Fast_H_sel=[]; % The selected channel matrix (empty array)
    
    % Initialize Alpha
    for j=1:Nr
        % H(j,:) = j-th row of channel matrix H. -> perspective of receiver
        % H(:,j) = j-th column of channel matrix H. -> pespective of transmit
        % alpha = hj*hj^H
        f=H(j,:);
        h=f';
        alpha=h'*h; % scalar
        Alpha=[Alpha alpha]; % record every alpha
    end
    
    % Cycle once to select an antenna
    for n=1:Lr
        % ----------------------- corresponding label -----------------------
        % Finding the j-th index means the index having the value that contributes the most to the entire channel matrix.
        [maxOfAlpha,index]=max(Alpha);
        % J_sel=[J_sel;J]; % The selected j-th transmit antenna index
        
        fullAntenna(index)=[]; % Remove the selected antenna
        Fast_H_sel=[Fast_H_sel;H(index,:)]; % The channel of the selected antenna
        % ----------------------- corresponding label -----------------------
        
        % Update B, and remove selected antenna
        if (n<Lr)
            f=H(index,:);
            h=f';
            alpha=Alpha(index);
            a=(B*h)/sqrt((Nt/SNR)+alpha);
            B=B-a*a';    
            
            % Remove the alpha corresponding to the selected antenna
            % L = {1, 2, ... Lt} - {J}
            % L = {L} - {J}
            % return L - {L} -> optimal maximum channel capacity index
            Alpha(index)=[];
            
            % The fullAntenna below has only the number of the optimal antenna subsets
            % alpha update
            for k=1:length(fullAntenna)
                Alpha(k)=Alpha(k)-(abs(a'*h))^2;
            end
        end
    end
    capacityOfSelected=log2(det(eye(Nt)+SNR/Nt*(Fast_H_sel'*Fast_H_sel)));
end