clc
clear all
tic

n = 4; %number of oscillators

nniter = 10; %number of niter monte-carlo simulations
niter = 100; %number of monte-carlo simulations

Vcn = 0; %initial conditions without perturbation
Ven = 0;
Iln = 0;

sigma = 5; %random number distribution magnitude

options = odeset('Events',@ClassEvent,'RelTol',1e-4,'AbsTol',1e-6);

lambdastep = 5; %how many steps you want lambda to run. 0 means only run for 
                %one desired lambda
data = zeros(lambdastep+1,8,nniter); %initializing data to be saved from MC sim

for k = 0:lambdastep %number of lambda steps you want to plot. 
    d = k/100; %lambda step size
    ilambda = -.49; %initial lambda

    lambda = ilambda + d;
    %p = [Vcn(0),Ven(0),ILn(0),Vcc,RL,L,C1,C2,Ree,Vee,beta,Ron,Vth,Vb,lambda,n]
    p = [5 200 21.8e-6 1e-9 1e-9 1000 -5 260 100 0.75 0 lambda n]; 

    transt = 8e-4; %transcient time
    for h = 1:nniter
        IPcount = 0; %intializing solution type count
        TWcount = 0;
        TW2count = 0;
        QPcount = 0;
        parfor kk = 1:niter %local pool parallel loop for montecarlo sim
            Vcn = 0; %initial conditions without perturbation
            Ven = 0;
            Iln = 0;

            sigma = 5; %random number distribution magnitude

            p0 = [Vcn; Ven; Iln]+(sigma*rand); %adds perturbation to IC
            p00 = [Vcn; Ven; Iln];

            for i = 1:n-1
                p0 = [p0;(p00+(sigma*rand))]; 
            end

            [t,Y] = ode23(@cc_oscillator2, [0:8e-8:transt], p0,[], p);

            inew = Y(end,:)'; %new initial conditions after transcient time

            [t1,Y1,te,ye,ie] = ode23(@cc_oscillator2, [0:8e-10:8e-6], inew, options, p);

%             figure()
%             for j = 1:n %plots all Vcn
%                 plot(t1,Y1(:,3*j-2),'LineWidth',1.5)
%                 hold on
%             end
%             hold off
%             title(['V_{Cn}: \lambda =' num2str(lambda)])
%             xlabel('t')
%             ylabel('V_{Cn}')
%     
%             grid on

%             figure() %3D plot
%     
%             plot3(Y1(:,1),Y1(:,4),Y1(:,7),'LineWidth',1.1)
%             xlabel('V_{C1}')
%             ylabel('V_{C2}')
%             zlabel('V_{C3}')
%             title(['Phase Plot V_{C1,2,3}: \lambda =' num2str(lambda)])
%         %     xlim([.2 5.3])
%         %     ylim([.2 5.3])
%         %     zlim([.2 5.3])
%             grid on

            for ii = 1 %uses event function to estimate periodicity 
                ind = find(ii==ie); 
                lmat = te(ind);
            end
            class1 = range(diff(sort(lmat)))*1e7; 
                        
            classIP = InPhase(Y1,n);
                        
            if mod(n,2) == 1
                if classIP < 1e-10 %if class close to 0, then IP.
                    IPcount = IPcount + 1;
                else
                    classTW = TravWave(t1,Y1,n);
                    if classTW >= 1 
                        if class1 < .1 %if periodicity is withing a certain 
                            %tolerance, we have single period solution travelling
                            %wave. Used to separate solutions with double
                            %periodicity.
                            TWcount = TWcount + 1;
                        else %if not, then it is quasi periodic or "other"
                            QPcount = QPcount + 1;
                        end
                    else 
                        QPcount = QPcount + 1;
                    end
                end
            else
                if classIP < 1e-10 %if class close to 0, then IP.
                    IPcount = IPcount + 1;
                else
                    classTW2 = TravWave2(Y1,n);
                    if classTW2 < 1e-3
                        TW2count = TW2count + 1;
                    else
                        classTW = TravWave(t1,Y1,n);
                        if classTW >= 1 
                            if class1 < .1 %if periodicity is withing a certain 
                                %tolerance, we have single period solution travelling
                                %wave. Used to separate solutions with double
                                %periodicity.
                                TWcount = TWcount + 1;
                            else %if not, then it is quasi periodic or "other"
                                QPcount = QPcount + 1;
                            end
                        else 
                            QPcount = QPcount + 1;
                        end
                    end
                end                
            end
        end
        
        lambda1 = lambda;
        lambda2 = lambda;
        lambda3 = lambda;
        lambda4 = lambda;
        IP = (IPcount/niter);
        TW = (TWcount/niter);
        TW2 = (TW2count/niter);
        QP = (QPcount/niter);

        pvalues = [IP TW TW2 QP lambda1 lambda2 lambda3 lambda4];
        data(k+1,:,h) = pvalues; 

        if IP == 0 %used for plotting to get rid of 0 values and make plot nicer
            IP = [];
            lambda1 = [];
        end

        if TW == 0 
            TW = [];
            lambda2 = [];
        end
        
        if TW2 == 0 
            TW2 = [];
            lambda3 = [];
        end

        if QP == 0
            QP = [];
            lambda4 = [];
        end

        %a dot means in phase. open circle is travelling wave. Asterisk is
        %other. an x is T/2 solution.

        figure(niter+1)
        plot(lambda1,IP,'k.',lambda2,TW,'ko',lambda3,TW2,'kx',lambda4,QP,'k*')
        hold on
        grid on
        ylabel('Probability')
        xlabel('Coupling Strength (\lambda)')
        title('CCO: Probability Density for Solution Types')
        ylim([0 1.1])
    end
    save('uniPDF_n4_2.mat','data')
end
% save('newPDF4.mat','data')
toc

function J = InPhase(Y,n)
    Vc = zeros(length(Y(end-10:end,1)),n);
    for ii = 1:n
        Vc(:,ii) = Y(end-10:end,3*ii-2);
    end
    
    Vcr = zeros(length(Vc(:,1)),n);
    
    for jj = 1:length(Vc(:,1))
        Vcr(jj) = range(Vc(jj,:));
    end

    J = range(Vcr); %this is used to classify as in phase or not. 
    %finds range between n Vcs. If close to 0, then it's in phase.
end

function I = TravWave(t1,Y,n)
    [~,peaks] = findpeaks(Y(:,1)); %finds local
    %maxima to estimate period 
    peaks = peaks(2:2:end-1);
    i_period = round(mean(diff(peaks))); 
%     t_period = mean(diff(t1(peaks)));
    i_shift = round(i_period/n);
    Tw_shift(:,1) = Y(end-50:end,1);


    for q = 2:n
        for qq = 0:(length(Tw_shift(:,1))-1)
            Tw_shift(qq+1,q) = Y((end-50+qq)-((q-1)*i_shift),1);
        end
    end

    Vcs = zeros(length(Y(end-50:end,1)),n);

    for tt = 1:n
        Vcs(:,tt) = Y(end-50:end,3*tt-2);
    end


    Tw_shift = Tw_shift(:,2:end);
    Vcsdiff = Vcs(:,2:end);

    ind = 1:1:(n-1);

    pind = perms(ind);

    p3 = zeros(length(Vcsdiff(:,1)),length(Vcsdiff(1,:)),length(pind(:,1)));
    check = zeros(length(Vcsdiff(:,1)),length(Vcsdiff(1,:)),length(pind(:,1)));
    rcheck = zeros(length(pind(1,:)),length(Vcsdiff(1,:)));
    rrcheck = zeros(length(pind(:,1)),1);

    for mm = 1:length(pind(:,1)) 
        for mmm = 1:n-1
            p3(:,mmm,mm) = Vcsdiff(:,pind(mm,mmm));
        end
        check(:,:,mm) = Tw_shift - p3(:,:,mm);
        for nn = 1:length(p3(1,:,1))
            rcheck(mm,nn) = sum(abs(check(:,nn,mm)));
        end
        rrcheck(mm) = min(rcheck(mm,:));
    end
    rrcheck = rrcheck';
    classTW = rrcheck < 1;
    I = sum(classTW);

end

function P = TravWave2(Y1,n)

    Vc = zeros(length(Y1(end-10:end,1)),n);
    for ii = 1:n
        Vc(:,ii) = Y1(end-10:end,3*ii-2);
    end

    VcTW2 = zeros(length(Y1(end-10:end,1)),2);

    VcTW2(:,1) = Vc(:,1);
    VcTW2(:,2) = Vc(:,3);

    Vcr2 = zeros(length(Vc(:,1)),n);
    for jj = 1:length(Vc(:,1))
        Vcr2(jj) = range(VcTW2(jj,:));
    end

    P = range(Vcr2);
end