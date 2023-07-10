function [A_contr,B_contr,C_contr,D_contr] = initVertPlasmaModel(LinearModel)
% Build model matrices
%
% x' = Ax+Bu+Ew'
% y  = Cx+Du+Fw
%
% The INPUT vector holds the voltages applied to the following cicuits:
% CS1, CS2U, CS2L, CS3U, CS3L , PF1-PF6, VS3 (in-vessel coils), VS1, Vpl
% Vpl is a voltage input on the plasma current equation
%
% The DISTUBANCES are:
% - betap*Ip
% - li*Ip
%
% The OUTPUT vector holds:
% - the currents in the CS/PF/VS3 coils (typically 12 elements)
% - the plasma current Ip
% - the gaps controlled by the old shape controller (6 gaps)
% - The gaps controlled by the XSC
% - the x-point descriptors
% - the plasma vertical position z_c
% 
% Here we are only interested in the channel from the in-vessel circuit 
% voltage to the in-vessel current and plasma vertical position outputs.

% PF indexes in the state and output vectors
PFidxOut  = LinearModel.PoloidalCircuits.OutputPosition(contains(LinearModel.PoloidalCircuits.Name,'IVS3'));

% PF indexes in the input vector
PFidxIn  = find(contains(LinearModel.InputsInfo.Name,'VS3'));

% Indexes for the other outputs 
ZCidxOut = LinearModel.OutputsInfo.OutputPosition(contains(LinearModel.OutputsInfo.Name,'CV-ZC'));

% Assemble state-space matrices
idxOut = [PFidxOut ZCidxOut];
idxIn  = PFidxIn;

A_contr = -inv(LinearModel.L)*LinearModel.R;
B_contr = LinearModel.L\LinearModel.S(:,idxIn);
C_contr = LinearModel.C(idxOut,:); 
D_contr = LinearModel.D(idxOut,idxIn);  


