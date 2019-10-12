clear all
data = importdata('RUN2unidir_colpitts_sigma.txt');

Niter = length(data(2,2:end));

N = 3:2:(2*Niter+1);

% max = length(data(2:end,1));

lam_data = data(2:end,1);

lamcounter = 1;
for jj = 2:length(lam_data)
    lamcounter = lamcounter + 1;
   if  (lam_data(1) == lam_data(jj))
       break
   end
end

max = lamcounter - 1;

lam_space = lam_data(1:max);

sigmaiter = length(lam_data)/max;

for kk = 0:sigmaiter-1
    for ii = 1:max 

        phasevec = data(1+(ii+kk*max),(2:end));
        xdata = log10(N);
        ydata = log10(phasevec);

        p = polyfit(xdata,ydata,1);
        yfit = xdata*p(1) + p(2);
        accuracy(ii) = norm(yfit - ydata);
        slope(ii) = p(1);

    end


    % xmaxa = max(N);
    % xmaxb = max(log10(N));


    xvec = lam_space;
    % xvec = 1:length(lam_space);
    figure(kk+2)
    plot(xvec,slope,'k-*')
    hold on
    plot(xvec,-.5*ones(size(lam_space)),'r--')
    xlabel('Coupling (\lambda)','fontsize',16)
    ylabel('Scaling Exponent (1/N^{m})','fontsize',16)
    title(['\sigma = ' num2str(data(1,kk+1)) ' (c++)'],'fontsize',16)
    
%     figure(1)

end