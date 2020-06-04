clear;
tic;
d = 50; % In km

lambda_l = 10; % Density of point process along the axes that gnnerates lines
lambda_c = 0.5; % Density of points points/km       (1 per km)
num_iterations = 100000;

% INITIALIZING VARIABLES
min_dist_d1(num_iterations) = 0;
min_dist_d2(num_iterations) = 0;
min_dist_x1(num_iterations) = 0;
min_dist_x2(num_iterations) = 0;

min_pathdist(num_iterations) = 0;


% Number of horizontal and vertical lines in each iteration
numhl = poissrnd(lambda_l*(2*d),1,num_iterations);
numvl = poissrnd(lambda_l*(2*d),1,num_iterations);


parfor j1=1:num_iterations
    Nhl = numhl(j1);
    Nvl = numvl(j1);
    rhoh1 = [0 -d+2*d*rand(1,Nhl) ];
    [rhoh2,indh2] = sort(abs(rhoh1));
    rhoh = rhoh1(indh2); % Generating horizontal lines
    rhov1 = [-d+2*d*rand(1,Nvl) ];%[ (rv+(d-rv)*urv).*sign(urv-.5) 0];
    [rhov2,indv2] = sort(abs(rhov1));% - rho1(end);
    rhov = rhov1(indv2);% Generating vertical lines
    
    
    len_chord = 2*d*ones(1, length(rhoh1)+length(rhov1)    );
    
    numc = poissrnd(lambda_c*len_chord);
    numc = numc + (numc==0);
    total_v = sum(numc);
    
    clens = cumsum(numc);
    idx=zeros(1,clens(end));idx2=idx; idx3=idx;
    idx([1 clens(1:end-1)+1]) = diff([0 len_chord]);
    len_c = cumsum(idx);
    idx2([1 clens(1:end-1)+1]) = diff([0 rhoh rhov]);
    rho_vec = cumsum(idx2);
    
    
    % GENERATING POINTS OF THE MPLCP
    x = -len_c/2 + len_c.*rand(1,total_v);
    numch = sum(numc(1:Nhl));
    numcv = sum(numc(Nhl + 1:end));
    nch = numc(1:Nhl);
    ncv = numc(Nhl+1:end);
    Ph = [ x(1:numch); rho_vec(1:numch) ].'; % Points on horizontal lines
    Pv = [ rho_vec(numch+1:total_v) ; x(numch+1:total_v) ].'; % Points on vertical lines
    
    P = [Ph; Pv]; % MPLCP points
    
    
    % COMPUTING THE SHORTEST PATH DISTANCE
    dr = min([abs(Ph(Ph(1:nch(1))>0))  d]);
    dl = min([abs(Ph(Ph(1:nch(1))<0))  d]);
    xr = min(abs(rhov(rhov>0)));
    xl = min(abs(rhov(rhov<0)));
    
    dlr = [dl dr];
    xlr = [xl xr];
    [min_dist_x1(j1), ind] = min(xlr);
    min_dist_x2(j1) = sum(xlr)-min_dist_x1(j1);
    min_dist_d1(j1) = dlr(ind);
    min_dist_d2(j1) = sum(dlr)-min_dist_d1(j1);
    
    wr = min( abs(P(:,1)-xr) + abs(P(:,2)) );
    wl = min( abs(P(:,1)+xl) + abs(P(:,2)) );
    
    min_pathdist(j1) = min([dr dl xr+wr xl+wl]);
    
end

%% THEORETICAL EXPRESSIONS
ll = lambda_l;
lc = lambda_c;

cdf_w1_1 = @(w) 1 - exp( -3*(lc+ll)*w +3*ll/2/lc*(1-exp(-2*lc*w)) );
cdf_w1_2 = @(w,x1,x2) 1 - exp( -3*(lc+ll)*w + ll/2/lc*(3 + 2*exp(-2*lc*x2) -exp(-2*lc*w) -4*exp(-lc*(x2+w)) ) );

cdf_w2_1 = @(w) 1 - exp( -3*(lc+ll)*w +3*ll/2/lc*(1-exp(-2*lc*w)) );
cdf_w2_2 = @(w,x1,x2) 1 - exp( -3*(lc+ll)*w + ll/2/lc*(3 + 2*exp(-2*lc*x1) -exp(-2*lc*w) -4*exp(-lc*(x1+w)) ) );

thy_cdf_Rm = @(rm) 1 - exp(-2*(ll+lc)*rm)...
    -2*ll*exp(-(ll+lc)*rm).*integral( @(x1)  (1-cdf_w1_1(rm-x1)).*exp(-(ll+lc)*x1), 0, rm, 'Arrayvalued', true)...
    -2*ll^2*integral( @(x2) exp(-(ll+lc)*x2).*(1-cdf_w2_1(rm-x2)).*integral( @(x1)  (1-cdf_w1_1(rm-x1)).*exp(-(ll+lc)*x1), ...
                rm-x2, x2, 'Arrayvalued', true), rm/2, rm, 'Arrayvalued', true)...
    -2*ll^2*integral( @(x2) exp(-(ll+lc)*x2).*integral( @(x1)  (1-cdf_w1_2(rm-x1,x1,x2)).*(1-cdf_w2_2(rm-x2,x1,x2)).*exp(-(ll+lc)*x1), ...
                0,rm-x2, 'Arrayvalued', true), rm/2, rm, 'Arrayvalued', true)...
    -2*ll^2*integral( @(x2) exp(-(ll+lc)*x2).*integral( @(x1)  (1-cdf_w1_2(rm-x1,x1,x2)).*(1-cdf_w2_2(rm-x2,x1,x2)).*exp(-(ll+lc)*x1), ...
                0,x2, 'Arrayvalued', true), 0, rm/2, 'Arrayvalued', true);

%% PLOT SIMULATION VS THEORY

resolution = 20;
rm_val = 0:max(min_pathdist)/resolution:max(min_pathdist);
thy_cdf = zeros(size(rm_val));
sim_cdf = zeros(size(rm_val));
for k1=1:length(rm_val)
    thy_cdf(k1) = thy_cdf_Rm(rm_val(k1)); % CDF from evaluating the theoretical expression
    sim_cdf(k1) = nnz(min_pathdist<=rm_val(k1))/length(min_pathdist); % CDF from Monte-Carlo simulations
end

figure; hold on
plot(rm_val, thy_cdf, '-r','linewidth',1.5);
plot(rm_val, sim_cdf, 'bo', 'linewidth', 1.5, 'MarkerFaceColor', [0 0 1]);
set(gca, 'FontSize',14)
xlabel('Distance $r_m$ (km)','FontSize',14,'Interpreter','latex');
ylabel('CDF $F_{R_m}(r_m)$','FontSize',14,'Interpreter','latex');
hg=legend('Theorem 2', 'Simulation', 'location','southeast');
set(hg, 'FontSize',14,'Interpreter','latex');
box on


toc