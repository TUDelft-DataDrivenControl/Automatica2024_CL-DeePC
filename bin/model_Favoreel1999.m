%% defining system from Favoreel 1999; SPC: Subspace Predictive Control
A = [ 4.40 1 0 0 0;
     -8.09 0 1 0 0;
      7.83 0 0 1 0;
     -4.00 0 0 0 1;
      0.86 0 0 0 0];
B = [0.00098 0.01299 0.01859 0.0033 -0.00002].';
C = eye(1,5);
D = 0;
K = [2.3 -6.64 7.515 -4.0146 0.86336].';

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

%% create plant models

% Kalman filter
plant = ss(A,[B K], C, [D eye(ny,nu)],[]);

