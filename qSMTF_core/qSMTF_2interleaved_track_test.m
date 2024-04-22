clear;
close all;
clc;

% the parameters for the simulated listener
% phi0 = [1.58, 2.81];
phi0 = [4, 1.2];
% phi0 = [1.2, 2.38];
% phi0 = [3.4, 0.71];

beta0 = 2;  % psychometric-function slope parameter 2
prob_lapse = 0.2;
lambda0 = 1-((1-prob_lapse) + 1/3*prob_lapse);  
pf0 = @( s, p, b,l ) 1/3 + (1-1/3-l) ./ ...
( 1 + exp( -b * ( s(:,1) - ( p(:,2) .* log( s(:,2) ./ p(:,1) ) ) ) ) );

qsmtf1 = qSMTF();
qsmtf2 = qSMTF();

ntrials = 40;
trial_seq = [ones(1,ceil(ntrials/2)), 2*ones(1,ceil(ntrials/2))];   % 1 for track1 and 2 for track2
trial_seq = trial_seq(randperm(ntrials));   % randomly interleave the two tracks

for itrial = 1:ntrials
    
    if trial_seq(itrial) == 1
        r = rand < pf0(qsmtf1.xnext, phi0, beta0, lambda0);  % r = 1 for correct, 0 for incorrect
        qsmtf1.update(r);

        figure(1);
        subplot(221);
        plot_x = qsmtf1.Rdepth_tot; %linspace(0, 35, 100);
        plot_y = qsmtf1.phi(end,2)*log(plot_x/qsmtf1.phi(end,1));
        plot_y0 = phi0(2)*log(plot_x/phi0(1));
        plot(plot_x, plot_y);
        hold on;
        plot(plot_x,plot_y0,'r--');
        scatter(qsmtf1.x(qsmtf1.r==1,2), qsmtf1.x(qsmtf1.r==1,1),'g');
        scatter(qsmtf1.x(qsmtf1.r==0,2), qsmtf1.x(qsmtf1.r==0,1),'r');
        plot(qsmtf1.xnext(2), qsmtf1.xnext(1), 'c+' );
        hold off;
        text(16, 10, num2str(qsmtf1.phi(end,:)));
        set(gca, 'XDir','reverse');
        xlabel('Modulation Depth (dB)');
        ylabel('Ripple Density (RPO)');
        title(['Trial Number: ' num2str(qsmtf1.n)]);
        axis([0 20 0 12]);
    
        subplot(222)
        a = linspace(qsmtf1.phi_lim(1,1), qsmtf1.phi_lim(2,1), qsmtf1.phi_grid(1)); % 
        b = linspace(qsmtf1.phi_lim(1,2), qsmtf1.phi_lim(2,2), qsmtf1.phi_grid(2)); % 
        surf(a, b, log(reshape(qsmtf1.phi_posterior, length(b), length(a))));
        xlabel('A');ylabel('B');
        shading interp;
        view(2);
        colorbar;
        clim([-75 -5]);
        axis([0 12 0 4]);
        hold on;
        plot(qsmtf1.phi(1:end,1), qsmtf1.phi(1:end,2), 'k.-' );
        plot(phi0(1), phi0(2), 'ko' );
        hold off;
        title( [ 'z-score = ' num2str(qsmtf1.zsc) ] );
        drawnow;

    elseif trial_seq(itrial) == 2
        r = rand < pf0(qsmtf2.xnext, phi0, beta0, lambda0);  % r = 1 for correct, 0 for incorrect
        qsmtf2.update(r);

        figure(1);
        subplot(223);
        plot_x = qsmtf2.Rdepth_tot; %linspace(0, 35, 100);
        plot_y = qsmtf2.phi(end,2)*log(plot_x/qsmtf2.phi(end,1));
        plot_y0 = phi0(2)*log(plot_x/phi0(1));
        plot(plot_x, plot_y);
        hold on;
        plot(plot_x,plot_y0,'r--');
        scatter(qsmtf2.x(qsmtf2.r==1,2), qsmtf2.x(qsmtf2.r==1,1),'g');
        scatter(qsmtf2.x(qsmtf2.r==0,2), qsmtf2.x(qsmtf2.r==0,1),'r');
        plot(qsmtf2.xnext(2), qsmtf2.xnext(1), 'c+' );
        hold off;
        text(16, 10, num2str(qsmtf2.phi(end,:)));
        set(gca, 'XDir','reverse');
        xlabel('Modulation Depth (dB)');
        ylabel('Ripple Density (RPO)');
        title(['Trial Number: ' num2str(qsmtf2.n)]);
        axis([0 20 0 12]);
    
        subplot(224)
        a = linspace(qsmtf2.phi_lim(1,1), qsmtf2.phi_lim(2,1), qsmtf2.phi_grid(1)); % 
        b = linspace(qsmtf2.phi_lim(1,2), qsmtf2.phi_lim(2,2), qsmtf2.phi_grid(2)); % 
        surf(a, b, log(reshape(qsmtf2.phi_posterior, length(b), length(a))));
        xlabel('A');ylabel('B');
        shading interp;
        view(2);
        colorbar;
        clim([-75 -5]);
        axis([0 12 0 4]);
        hold on;
        plot(qsmtf2.phi(1:end,1), qsmtf2.phi(1:end,2), 'k.-' );
        plot(phi0(1), phi0(2), 'ko' );
        hold off;
        title( [ 'z-score = ' num2str(qsmtf2.zsc) ] );
        drawnow;
    end
end

phi_final = (qsmtf1.phi(end,:) + qsmtf2.phi(end,:))/2;
plot_y_final = phi_final(2)*log(plot_x/phi_final(1));
figure(1)
subplot(221)
hold on;
plot(plot_x, plot_y_final, 'm');
hold off;
subplot(223)
hold on;
plot(plot_x, plot_y_final, 'm');
hold off;
