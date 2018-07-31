function [value, isterminal, direction]=ClassEvent(t,y,p)
A = 2;
n = p(13);

% Vc = zeros(1,n);
% for i = 1:n
% Vc(i) = y(3*n-2);
% end

value = [A-y(1:3:(3*n))];
% direction  = 1;
% isterminal = 0;


direction  = ones(n,1);
isterminal = zeros(n,1);
% direction  = 1;
% isterminal = 0;


end