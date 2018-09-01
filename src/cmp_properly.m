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
S = 300;
T_primal = [];  % cvx_cpu_times
T_dual = [];
for s = 2:S     % s - number discrete states of the input random variables
    n = s^2;    % n - number of states of the joint probability

    P_x = rand(1,s);    P_x = P_x/sum(P_x);
    P_y = rand(1,s);    P_y = P_y/sum(P_y);
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
        variable x(n)
        cost = D' * x;

        minimize (cost)
        subject to 
            A*x == b
              x >= 0
    cvx_end
    T_primal = [T_primal cvx_cputime];
    
    
    cvx_begin
        variable x(2*s)
        cost = b' * x;

        maximize (cost)
        subject to 
            A'*x <= D
    cvx_end
    T_dual = [T_dual cvx_cputime];
end
dlmwrite('primal_times_300.txt',T_primal)
dlmwrite('dual_times_300.txt',T_dual)


figure
X = linspace(2,S,S-1);
plot(X,T_primal,'b-*')

hold on
grid on
xlabel('number of states of the input variables')
ylabel('time')
plot(X,T_dual, 'r-*')  
legend('primal','dual')
legend('show')
