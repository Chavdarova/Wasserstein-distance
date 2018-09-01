% Copyright (c) 2008 Idiap Research Institute, http://www.idiap.ch/
% Written by Tatjana Chavdarova <tatjana.chavdarova@idiap.ch>

% This file is part of Wasserstein Distance project.

% Wasserstein Distance is free software: you can redistribute it and/
% or modify it under the terms of the GNU General Public License 
% version 3 as published by the Free Software Foundation.

% Wasserstein Distance is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this project. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 5;      % number discrete states of the input random variables
n = s^2;    % number of states of the joint probability (and the 
            % dimension of our optimization variable)


%% Example 1 in table 1.
% P_x = [2 2 2 2 2];    P_x = P_x/sum(P_x)  
% P_y = [2 2 2 2 2];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +4.45945e-13

%% Example 2 in table 1.
% P_x = [10 8 6 4 2];    P_x = P_x/sum(P_x)
% P_y = [2 4 6 8 10];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +1.33333

%% Example 3 in table 1.
% P_x = [10 1 1 1 1];    P_x = P_x/sum(P_x)
% P_y = [1 1 1 1 10];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +2.57143

%% Example 4 in table 1.
% P_x = [20 5 1 30 1];    P_x = P_x/sum(P_x)
% P_y = [1 5 10 0 10];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +1.04656

%% Example 5 in table 1.
% P_x = [0 0 1 0 1];    P_x = P_x/sum(P_x)
% P_y = [1 1 0 1 0];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +1.66667

%% Example 6 in table 1.
P_x = [0 0 1 0 0];    P_x = P_x/sum(P_x)
P_y = [1 0 0 0 0];    P_y = P_y/sum(P_y)
% Optimal value (cvx_optval): +2

b = [P_x P_y]';

A_x = [];
for i = 1:s
    B = zeros(s);
    B(i, :) = 1;
  A_x = cat(2, A_x, B);
end
A_y = repmat(diag(ones(s,1)), 1, s);
A = [A_x; A_y];

D = zeros(s);
for i = 1:s
    for j = 1:s
        D(i,j) = abs(i-j);
    end
end
D = reshape(D, n,1);

cvx_begin
    variable x(2*s)
    cost = b' * x;
    
    maximize (cost)
    subject to 
        A'*x <= D
cvx_end