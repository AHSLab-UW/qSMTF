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

qsmtf = qSMTF();

ntrials = 40;
for itrial = 1:ntrials
    if qsmtf.repeat_flag == 1
        disp('stimulus repeated.');
    end
    r = rand < pf0(qsmtf.xnext, phi0, beta0, lambda0);  % r = 1 for correct, 0 for incorrect
    qsmtf.update(r);

    figure(1);
    subplot(121);
    plot_x = qsmtf.Rdepth_tot; %linspace(0, 35, 100);
    plot_y = qsmtf.phi(itrial+1,2)*log(plot_x/qsmtf.phi(itrial+1,1));
    plot_y0 = phi0(2)*log(plot_x/phi0(1));
    plot(plot_x, plot_y);
    hold on;
    plot(plot_x,plot_y0,'r--');
    scatter(qsmtf.x(qsmtf.r==1,2), qsmtf.x(qsmtf.r==1,1),'g');
    scatter(qsmtf.x(qsmtf.r==0,2), qsmtf.x(qsmtf.r==0,1),'r');
    plot(qsmtf.xnext(2), qsmtf.xnext(1), 'c+' );
    hold off;
    text(16, 10, num2str(qsmtf.phi(end,:)));
    set(gca, 'XDir','reverse');
    xlabel('Modulation Depth (dB)');
    ylabel('Ripple Density (RPO)');
    title(['Trial Number: ' num2str(itrial)]);
    axis([0 20 0 12]);

    subplot(122)
    a = linspace(qsmtf.phi_lim(1,1), qsmtf.phi_lim(2,1), qsmtf.phi_grid(1)); % 
    b = linspace(qsmtf.phi_lim(1,2), qsmtf.phi_lim(2,2), qsmtf.phi_grid(2)); % 
    surf(a, b, log(reshape(qsmtf.phi_posterior, length(b), length(a))));
    xlabel('A');ylabel('B');
    shading interp;
    view(2);
    colorbar;
    clim([-75 -5]);
    axis([0 12 0 4]);
    hold on;
    plot(qsmtf.phi(1:itrial+1,1), qsmtf.phi(1:itrial+1,2), 'k.-' );
    plot(phi0(1), phi0(2), 'ko' );
    hold off;
    title( [ 'z-score = ' num2str(qsmtf.zsc) ] );
    drawnow;


end