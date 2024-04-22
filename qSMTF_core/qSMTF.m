classdef qSMTF <handle
    
    properties
        %basic properties
        x    % stimulus parameter
        xnext % the next stimulus based on the previous data;
        Rdepth_tot = [0.5:0.5:2, 3:1:20]';      % spectral modulation depths
        RPO_tot = [ 0.125 0.149 0.177 0.210 0.250 0.297 ...     % spectral modulation rate
	        0.354 0.420 0.500 0.595 0.707 0.841, ...
	        1 1.189 1.414 1.682 2 2.378 2.828 3.364 ...
	        4 4.757 5.657 6.727 8 9.514 11.314 13.454 16 19.027 20]';
        stim_par    % all possible stimuli
        
        
        % parameter space
        phi_lim = [0.2 0.1; 12 4];
        phi_grid = [100 120];
        A
        B
        phi_prior
        phi_posterior
        phi % model parameters
        beta = 5;
        lambda = 0.1;
        nafc = 3;
        pf

        % response and control parameters
        r   %the listener's responses (correctness) 
        n=0;   %the trial number
        repeat_flag = 0;
        zsc
        
    end
    
    methods

        %constructor
        function qsmtf = qSMTF()
            reset(qsmtf);
        end  
        
        % INITIALIZATION
        function reset(qsmtf)
            qsmtf.x;
            qsmtf.r;
            qsmtf.n = 0;

            % psychometric function
            qsmtf.pf = @( s, p, b, l) 1/qsmtf.nafc + (1-1/qsmtf.nafc-l) ./ ...
	        ( 1 + exp( -b * ( s(:,1) - ( p(:,2) .* log( s(:,2) ./ p(:,1) ) ) ) ) );

            % stimulus space
            [stim_x, stim_y] = meshgrid(qsmtf.RPO_tot, qsmtf.Rdepth_tot);
            qsmtf.stim_par = [stim_x(:), stim_y(:)]; % all possible stim combinations
            
            % parameter space
            a = linspace(qsmtf.phi_lim(1,1), qsmtf.phi_lim(2,1), qsmtf.phi_grid(1)); % 
            b = linspace(qsmtf.phi_lim(1,2), qsmtf.phi_lim(2,2), qsmtf.phi_grid(2)); % 
            [A1, B1] = meshgrid(a,b);
            qsmtf.A = A1(:);
            qsmtf.B = B1(:); % s.t. [A,B] is all the possible phi combinations,

            % prior distributions
            % s.t. pf( ..., [A,B], ...) is likelihood over all possible phi
            phi_prior_init = ones(size(qsmtf.A))/length(qsmtf.A); % initial prior is uniform over phi
            qsmtf.phi_prior = phi_prior_init;
            qsmtf.phi(1,:) = [sum(qsmtf.A.*qsmtf.phi_prior) sum(qsmtf.B.*qsmtf.phi_prior)];

            % find the stimulus for the first trial
            findx(qsmtf);
        end
        
        % ITERATION
        
        %Update posterior and xnext based on the previous posterior and the
        %new response r.
        function update(qsmtf, r , x)
            qsmtf.n = qsmtf.n + 1;
            qsmtf.r(qsmtf.n,:) = r;
            if nargin == 2
                qsmtf.x(qsmtf.n,:) = qsmtf.xnext;
            elseif nargin == 3
                qsmtf.x(qsmtf.n,:) = x;
            end
            
            % update model parameters
            % calculate the likelihood of the response
            p_r = qsmtf.pf(qsmtf.x(qsmtf.n,:), qsmtf.phi(qsmtf.n,:), qsmtf.beta, qsmtf.lambda);
            r_likelihood = p_r^(r)*(1-p_r)^(1-r);
    
            % if r is not expected from the model, we will repeat collecting
            % data at this place
            % step 3: re-fit the model
            if qsmtf.n > 10 && r_likelihood <1/qsmtf.nafc
                qsmtf.repeat_flag = 1;
            else
                qsmtf.repeat_flag = 0;
            end
            
            likelihood = qsmtf.pf(qsmtf.x(qsmtf.n,:), [qsmtf.A, qsmtf.B], qsmtf.beta, qsmtf.lambda);
            if r ~= 1, likelihood = 1-likelihood; end
            qsmtf.phi_posterior = qsmtf.phi_prior .* likelihood; % update the heatmap
            qsmtf.phi_posterior = qsmtf.phi_posterior ./ sum(qsmtf.phi_posterior);   % normalize the probability distribution
        
            % find the model parameters with the maximum likelihood
            % phi(itrial + 1, :) = [sum(A.*phi_posterior) sum(B.*phi_posterior)]; 
            [~, maxL_idx] = max(qsmtf.phi_posterior);
            qsmtf.phi(qsmtf.n + 1, :) = [qsmtf.A(maxL_idx) qsmtf.B(maxL_idx)];     % new fitted model parameters
            qsmtf.phi_prior = qsmtf.phi_posterior;  % the updated posterior distribution becomes the prior distribution for the next trial
                
            %Find the next signal strength
            findx(qsmtf);
        end

        %Stimulus selection
        function findx(qsmtf)
            if qsmtf.repeat_flag == 1
                qsmtf.xnext = qsmtf.x(end,:);     % repeat the last stimulus
            else
                % stimulus selection based on minimizing entropy
	            I = zeros(length(qsmtf.stim_par),1); % info (negentropy) difference between yes and no
	            for istim = 1:length(qsmtf.stim_par)
		            p = qsmtf.pf( qsmtf.stim_par(istim,:), [qsmtf.A,qsmtf.B], qsmtf.beta, qsmtf.lambda); % prob of all phi given this stim
		            p0 = qsmtf.phi_prior .* (1-p);
		            p0 = p0./sum(p0); % conditional for no
		            p1 = qsmtf.phi_prior .* p;
		            p1 = p1./sum(p1); % conditional for yes
		            p_pred = qsmtf.pf( qsmtf.stim_par(istim,:), qsmtf.phi(qsmtf.n + 1, :), qsmtf.beta, qsmtf.lambda); % prob of current phi given current stim
		            I(istim) = -(1-p_pred) * sum( p0.*log(p0) ) - p_pred * sum( p1.*log(p1) ); % info(no) - info(yes)
	            end
            % 	disp((min(I)-mean(I))/std(I));
	            qsmtf.zsc = (min(I)-mean(I))/std(I);  % candidate metric for determine outlier?
            
                % if the model is already reasonably converged, then choose the
                % stimulus that will lead to the minimum expected entropy (maximum
                % expected information gain). Otherwise, we will select a random
                % modulation depth, and then sample the RPO with the minimum expected
                % entropy for that modulation depth.
            
	            if qsmtf.zsc < -3 % z-score of min is large, so...
		            % minimize w.r.t entire stimulus space
                    % disp('zsc < -3');
		            [~, min_entropy_idx] = min(I);
		            qsmtf.xnext = qsmtf.stim_par(min_entropy_idx,:);
	            else
		            % otherwise, perturb by minimizing w.r.t randomly selected modulation depth
		            depth_idx = find(qsmtf.stim_par(:,2) == qsmtf.Rdepth_tot(randi(length(qsmtf.Rdepth_tot))));
		            [~, min_entropy_idx] = min(I(depth_idx));
		            qsmtf.xnext = qsmtf.stim_par(depth_idx(min_entropy_idx),:);
                end
            end
        end

    end

end

%eof