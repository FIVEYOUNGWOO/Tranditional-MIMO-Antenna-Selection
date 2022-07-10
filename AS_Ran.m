% Random-based antenna selection in MIMO System

function  [capacityOfSelected, Rand_H_sel]=AS_Ran(Nr,Nt,Lr,SNR,H,fullAntenna)
if(Lr==Nr)
    % Selected channel capacity
    capacityOfSelected=log2(det(eye(Nt)+SNR/Nt*(H'*H))) ;
else
    
    % Selected channel matrix
    Rand_H_sel=[];
    
    % Cycle to select an antenna
    for n=1:Lr
        
        % for select antenna elements
        randomIndex=randi([1, length(fullAntenna)],1,1);
        
        % The channel of the selected antenna
        Rand_H_sel=[Rand_H_sel;H(randomIndex,:)];
        
        % Remove the antenna that has been selected
        fullAntenna(randomIndex)=[];
    end
    
    % Selected channel capacity
    capacityOfSelected=log2(det(eye(Nt)+SNR/Nt*(Rand_H_sel'*Rand_H_sel))) ;
end