# qSMTF
A Matlab toolbox for running SMTF experiments

The basic structure of building an adaptive track is:
```
qsmtf = qSMTF();
for itrial = 1:ntrials
    % getting a response r from the participant
    qsmtf.update(r);
end
phi = qsmtf.phi(end,:);
```
 
The key variables are:
```
qsmtf.x  % a list of all stimuli [modulation rate (RPO), modulation depth (dB)]
qsmtf.r  % a list of all responses
qsmtf.phi             % a list of all estimated model parameters (intercept and slope)
qsmtf.xnext         % the stimulus for the next trial, selected by the algorithm
qsmtf.posterior  % the posterior parameter distribution
```

The algorithm is currently implemented to automatically repeat a stimulus when a questionable response has been received.

I’ve provided two test scripts. One is the unit test for the qSMTF. The other provides an example of two randomly interleaved adaptive tracks. For these test experiments, I’ve set the lapse rate to 20% (20% of trials with random button presses). 
