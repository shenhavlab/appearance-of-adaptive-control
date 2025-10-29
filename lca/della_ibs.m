%% Fit LCA
function della_ibs(ibs_maxtime)
% clear;

addpath(genpath('./bads'))
addpath(genpath('./ibs'))

n_inits = 40;

max_bads_iter = 100
MaxFunEvals = 1000

ibs_maxtime = double(ibs_maxtime)


global_disc = 0
n_discs = 10
dt = .005

center_rt = true


savename = sprintf('times-%d_centerRT_ibs-2',ibs_maxtime*1000)

%% prep data
PATH = 'DATA_PATH';
D=readtable('speeded_behavioral_data.csv');
D = D(~D.responseID==0, :); % removing 3 trials where there was no response
D.responseID(D.responseID==-1) = 0;
tablesize = height(D);

pts = unique(D.subjectID);
n_pts = length(pts)

%% convert to accuracy code
Vs = [ D.ExpectedValue_right, D.ExpectedValue_left];
[~,max_idx] = max(Vs,[],2);
[~,min_idx] = min(Vs,[],2);

D.acc = (D.responseID+1) == max_idx;
D.rt = D.ResponseTime_seconds;
D.maxRT = zeros(size(D.rt)) + .750;
D.acc(isnan(D.rt))=0;% treat misses as errors


if center_rt
    % center = @(x) x-nanmean(x)
    grand_mean = nanmean(D.rt);

    for pp = 1:n_pts

        mRT = nanmean(D.rt(D.subjectID==pts(pp)));

        D.rt(D.subjectID==pts(pp)) = D.rt(D.subjectID==pts(pp)) - mRT;
        D.maxRT(D.subjectID==pts(pp)) = D.maxRT(D.subjectID==pts(pp)) - mRT;

    end

    D.rt = D.rt + grand_mean;
    D.maxRT = D.maxRT + grand_mean;

end




assert(mean(D.acc) > .5)

[D.CorrMaxVal, D.CorrMinVal, D.ErrMaxVal, D.ErrMinVal,D.CorrVal, D.ErrVal] = deal(D.ExpectedValue_right);
for tt = 1:tablesize
    if max_idx(tt) == 1

        D.CorrVal(tt) = D.ExpectedValue_right(tt);
        D.ErrVal(tt) = D.ExpectedValue_left(tt);

        D.CorrMaxVal(tt) = D.MaxValue_right(tt);
        D.CorrMinVal(tt) = D.MinValue_right(tt);
        D.ErrMaxVal(tt) = D.MaxValue_left(tt);
        D.ErrMinVal(tt) = D.MinValue_left(tt);

    else

        D.CorrVal(tt) = D.ExpectedValue_left(tt);
        D.ErrVal(tt) = D.ExpectedValue_right(tt);

        D.CorrMaxVal(tt) = D.MaxValue_left(tt);
        D.CorrMinVal(tt) = D.MinValue_left(tt);
        D.ErrMaxVal(tt) = D.MaxValue_right(tt);
        D.ErrMinVal(tt) = D.MinValue_right(tt);

    end

end

assert(all(D.CorrVal >= D.ErrVal))


%% Get VD and OV

D.OV = D.ExpectedValue_right+D.ExpectedValue_left;
D.discOV = D.OV > median(D.OV);

D.VD = abs(D.ExpectedValue_right-D.ExpectedValue_left);
D.discVD = D.VD > median(D.VD);


%% compute quantiles within each OV condition
cc =1;
quants = nan(4,n_discs+1);
mean_acc = nan(4,1);

[grp_OV, grp_VD, grp_maxRT, grp_RT, grp_acc, grp_CorrMaxVal, grp_CorrMinVal, grp_ErrMaxVal,grp_ErrMinVal] = deal([]);

for ov = 0:1
    for ac = 0:1

        Dg = D(D.discOV==ov & D.acc==ac,:);

        % create dmat
        grp_OV          = [grp_OV; Dg.discOV];
        grp_VD          = [grp_VD; Dg.discVD];
        grp_acc         = [grp_acc; Dg.acc];

        mean_acc(cc)    = nanmean(Dg.acc);

        quants(cc,:)    = quantile(Dg.rt, linspace(0,1,n_discs+1));
        % quants(cc,:)    = linspace(nanmin(D.rt), nanmax(D.rt),n_discs+1);
        grp_RT          = [grp_RT; discretize(Dg.rt, quants(cc,:))];

        grp_CorrMaxVal  = [grp_CorrMaxVal; Dg.CorrMaxVal];
        grp_CorrMinVal  = [grp_CorrMinVal; Dg.CorrMinVal];
        grp_ErrMaxVal   = [grp_ErrMaxVal; Dg.ErrMaxVal];
        grp_ErrMinVal   = [grp_ErrMinVal; Dg.ErrMinVal];

        grp_maxRT       = [grp_maxRT; Dg.maxRT];


        cc=cc+1;
    end
end

% set missed responses to their own bin
grp_RT(isnan(grp_RT)) = 0;


if global_disc

    disp('RUNNING GLOBAL')
    quants = quantile(D.rt, linspace(0,1,n_discs+1));
    grp_RT = discretize(D.rt, quants);
    grp_RT(isnan(grp_RT)) = 0;

end

quants


%% convert to matrix

designMat = [...
    grp_OV, grp_acc,...
    grp_CorrMaxVal/80, grp_CorrMinVal/80, grp_ErrMaxVal/80, grp_ErrMinVal/80, ...
    grp_VD, grp_maxRT...
    ];

respMat = grp_RT; 


%% LCA parameters

badsopt = bads('defaults');
badsopt.UncertaintyHandling = true;
badsopt.SpecifyTargetNoise = true;
badsopt.NoiseFinalSamples   = 30;
badsopt.MaxIter = max_bads_iter;
badsopt.Plot = 'off';
badsopt.MaxFunEvals=MaxFunEvals;
badsopt

param_names = {'t0', 'vin', 'vratio', 'leak', 'inhib', 'bound', 'collapse', 'sigma'};
x0 = [          0,     6,    .75,      1.75,   .03,     3,       1,          .8];
LB = [          0,     1,    .01,      0,      0,       2,       0,          .1];
UB = [          .300,  30,   .99,      5,      5,       20,      2,          10];
PLB = [         .200,  1,    .6,       .5,     .5,      3,       .5,         0.5];
PUB = [         .255,  10,   .9,       2,      2.0,     10,      1.5,        1.5];

n_param = length(param_names)

%% IBS options

ibsopt = ibslike('defaults');
ibsopt.Nreps = 10;
ibsopt.ReturnStd = true;
ibsopt.MaxTime = ibs_maxtime;
% ibsopt.MaxIter=5;
ibsopt

%% fit options

llopt = struct;
llopt.quants = quants;
llopt.lates = mean(grp_RT==0);
llopt.accs = mean_acc;
llopt.global_disc = global_disc;
llopt.maxRT = nanmax(D.rt);
llopt.dt = dt;
llopt

%% create functions

inner_fun = @(x,dmat) lca_sim_ibs(x, dmat, llopt)
outer_fun = @(x) ibslike(inner_fun,x,respMat, designMat, ibsopt)


%% test inner


[resp,rt] = lca_sim_ibs(x0, designMat, llopt);

n_tests = 200
tic
for ii = 1:n_tests
    ll = inner_fun(x0, designMat);
end
inner_dur = toc;
fprintf('\ninner dur = %.2g s/iter\n', inner_dur/n_tests)


%% test outer


[NLOGL,NLOGLVAR,EXITFLAG,OUTPUT] = outer_fun(x0)



%% run fit

% [x,fval,exitflag,output,optimState,gpstruct] = bads(outer_fun,x0,LB,UB,PLB,PUB,[],badsopt)

rep_fval = nan(n_inits,1);
rep_x = nan(n_inits, n_param);
[rep_exitflag, rep_output, rep_optimState, rep_gpstruct] = deal(cell(n_inits,1));


time_on = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'))
for rr = 2:n_inits
    tic

    if rr ==1
        [rep_x(rr,:),rep_fval(rr),rep_exitflag{rr},rep_output{rr},rep_optimState{rr},rep_gpstruct{rr}] = bads(outer_fun,x0,LB,UB,PLB,PUB,[],badsopt);
    else
        rep_x0 = unifrnd(PLB, PUB);
        [rep_x(rr,:),rep_fval(rr),rep_exitflag{rr},rep_output{rr},rep_optimState{rr},rep_gpstruct{rr}] = bads(outer_fun,rep_x0,LB,UB,PLB,PUB,[],badsopt);
    end

    % print
    tfit=toc;
    fprintf("\niter=%d // fval=%.2g // dur=%.3g min", rr, rep_fval(rr), tfit/60)

    % best
    [best_fval,best_fit]    = min(rep_fval)
    best_x          = rep_x(best_fit,:)
    best_exitflag   = rep_exitflag{best_fit};
    best_output     = rep_output{best_fit};
    best_optimState = rep_optimState{best_fit};
    best_gpstruct   = rep_gpstruct{best_fit};

    % save
    save(sprintf('fit-results/%s_%s', time_on, savename))

end

best_x
best_exitflag
best_output
best_optimState
best_gpstruct

