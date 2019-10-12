clear all

% fileID2 = fopen('lambdaOUT.bin');
% lam_space = fread(fileID2,'double');
% 
% max = length(lam_space);
% slope_space = zeros(1,max);
max =1;


for ii = 1:max
%     fileID = fopen(['phaseOUT_iter' num2str(ii-1) '.bin']);
    fileID = fopen(['phaseOUT_iter' num2str(100) '.bin']);
    phase_space = fread(fileID,'double');


    Nsize = length(phase_space);

    N = 3:2:(2*Nsize+1);

    xdata = log(N);
    ydata = log(phase_space');

    p = polyfit(xdata,ydata,1)

    slope_space(ii) = p(1);
    figure()
    plot(xdata,ydata,'b*')
    
    fileID = [];
    phase_space = [];
%     
end
fs = 16;
figure()
plot(lam_space, slope_space, 'k-*')
hold on
plot(lam_space,-.5*ones(size(lam_space)),'r')

xlabel('Coupling (\lambda)','fontsize',16)
ylabel('Scaling Exponent (1/N^{m})','fontsize',16)