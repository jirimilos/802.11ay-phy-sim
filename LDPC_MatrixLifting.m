function [ayLDPCLiftedMatrix] = LDPC_MatrixLifting(adLDPCMatrix, LiftingMatrix)
%LDPC_MatrixLifitng Function provides 2nd lifitng of the original LDPC matrix
% Author:   Jiri Milos, DREL FEEC BUT, 2020

% test:
% LiftingMatrix = [1 0 1; 0 0 0; 1 0 1];
% adLDPCMatrix = [29 30 0; 37 31 18; 25 22 4];

ayLDPCLiftedMatrix = -1*ones(2*size(adLDPCMatrix));
bMatrix = [1 -1; -1 1];
for ir = 1:size(adLDPCMatrix,1)
    adLDPCMatrix_row = adLDPCMatrix(ir, :);
    LiftingMatrix_row = LiftingMatrix(ir,:);
    for ic = 1:size(adLDPCMatrix,2)
        adLDPCMatrix_item = adLDPCMatrix_row(ic);
        LiftingMatrix_item = LiftingMatrix_row(ic);
        if LiftingMatrix_item == -1 % if Lifitng item equals -1 then leave all corresponding items equaling to -1
            continue;
        else % otherwise
            if adLDPCMatrix_item == 0
                ayLDPCLiftedMatrix(((ir-1)*2)+[1 2],((ic-1)*2)+[1 2]) = [-1 0; 0 -1];
            else
                ayLDPCLiftedMatrix(((ir-1)*2)+[1 2],((ic-1)*2)+[1 2]) = adLDPCMatrix_item*circshift(bMatrix, LiftingMatrix_item);
            end
        end
    end
end

% set all-zero submatrices to -1 (NaN)
% ayLDPCLiftedMatrix(ayLDPCLiftedMatrix < 0) = -1;
ayLDPCLiftedMatrix(ayLDPCLiftedMatrix < 0) = NaN;

end

