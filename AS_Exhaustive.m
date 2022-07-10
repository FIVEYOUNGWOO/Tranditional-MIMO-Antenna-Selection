function  [capacityOfSubsetMax, Ex_H_sel]=AS_Exhaustive(Nr,Nt,Lr,SNR,H,antennaSubset)

    % Record the maximum channel capacity of all subset when the number of antennas selected is Lr
    capacityOfSubsetMax=0;
    Ex_H_sel = [];
    
    % The number of antennas selected is Lr.
    % The loop produces one subset of channels at a time.
    % With a total of nchoosek (Nr,Lr) subset.
    for k=1:nchoosek(Nr,Lr) 
        % Subset
        indexOfChannel=antennaSubset(k,:);
        
        % Selected channel matrix H
        Ex_H_sel=H(indexOfChannel,:);
        
        % Capacity of subset
        capacityOfSubset=log2(det(eye(Nt)+SNR/Nt*(Ex_H_sel'*Ex_H_sel)));
        
        % Maximum channel capacity among all subsets
        if(capacityOfSubset>capacityOfSubsetMax)
            capacityOfSubsetMax=capacityOfSubset;
        end
    end
