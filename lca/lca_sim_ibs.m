%% LCA simulation function
function [resps, rts, choices] = lca_sim_ibs(x, dmat, llopt)


% param_names = { 't0', 'vin', 'vratio', 'leak', 'inhib', 'bound', 'collapse', 'sigma'};
%
% designMat = [...
%     grp_OV, grp_acc,...
%     grp_CorrMaxVal, grp_CorrMinVal, grp_ErrMaxVal, grp_ErrMinVal...
%     ];
%
% respMat = [...
% grp_RT, grp_acc
%     ];

%% handle parameters


t0       = x(1);
vin      = x(2);
vratio   = x(3);
leak     = x(4);
inhib    = x(5);
bound    = x(6);
collapse = x(7);
sigma    = x(8);


%% constants

dt  = llopt.dt;
sdt = sigma*sqrt(dt);

n_accumulators  = 2;
n_trials        = size(dmat,1);
steps           = t0:dt:llopt.maxRT;
nstep           = length(steps);
n_active        = n_trials;

% preallocate
x           = zeros(n_accumulators, n_trials);
active      = true(1, n_trials);
active_mask = true(n_accumulators,n_trials);
choices     = nan(n_trials,1);
rts         = nan(n_trials,1);
maxRT       = dmat(:,8);



% process inputs
%'CorrMaxVal', 'CorrMinVal', 'ErrMaxVal', 'ErrMinVal'}
W_in = dt*vin*[vratio, (1-vratio), 0, 0; 0, 0, vratio, (1-vratio)];
drive_in = W_in*dmat(:,3:6)';


% bounds
bound_t =  bound - bound*collapse*steps;

% inhibition
W_x = eye(2) + dt*[-leak, -inhib; -inhib, -leak];


for t = 1:nstep

    if ~n_active, break; end  % Early termination saves RNG calls

    % LCA dynamics (vectorized)
    x(:,active) = max(0,...
        W_x*x(:,active) +...
        drive_in(:,active) +...
        randn(n_accumulators, n_active) * sdt...
        );

    % Check thresholds
    hit_trial = any((x >= bound_t(t)) & active_mask,1);
    if any(hit_trial)
        [~,choices(hit_trial)] = max(x(:,hit_trial),[],1);
        rts(hit_trial) = steps(t);
        active(hit_trial) = false;
        active_mask(:,hit_trial) = false;
        n_active = sum(active);
    end


end



%% compute quantiles

if llopt.global_disc

    disc_rts = discretize(rts, llopt.quants);
    disc_rts(active) = 0; % missed resposnes
    disc_rts(isnan(disc_rts)) = -1; % early responses

    resps = [disc_rts, choices];

else % discretize by OV condition

    disc_rts = rts;
    cc=1;
    for ov = 0:1
        for ac = 0:1

            disc_rts(dmat(:,1)==ov & dmat(:,2)==ac) = discretize(disc_rts(dmat(:,1)==ov & dmat(:,2)==ac), llopt.quants(cc,:));
            cc=cc+1;

        end
    end

    % missed responses
    disc_rts(active' | (rts > maxRT)) = 0; 
    % choices(disc_rts==0) = 0;

    % early responses
    disc_rts(isnan(disc_rts)) = -1; 

    % combine
    % resps = [disc_rts, choices];
    resps = disc_rts;

end

end

