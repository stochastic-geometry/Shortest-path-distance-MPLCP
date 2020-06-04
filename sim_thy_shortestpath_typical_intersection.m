clear;
tic;
d = 50; % In km
 
lambda_l = 10; % Density of point process along the axes that gnnerates lines
lambda_c = 0.5; % Density of points points/km       (1 per km)
num_iterations = 100000;

min_pathdist(num_iterations) = 0;

numhl = poissrnd(lambda_l*(2*d),1,num_iterations); % Number of horizontal lines in each iteration
numvl = poissrnd(lambda_l*(2*d),1,num_iterations); % Number of vertical lines in each iteration

for j1=1:num_iterations
    Nhl = numhl(j1);
    Nvl = numvl(j1);
    rhoh1 = [0 -d+2*d*rand(1,Nhl) ];
    [rhoh2,indh2] = sort(abs(rhoh1));
    rhoh = rhoh1(indh2); % Generating horizontal lines
    rhov1 = [0 -d+2*d*rand(1,Nvl) ];%[ (rv+(d-rv)*urv).*sign(urv-.5) 0];
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

    P = [Ph; Pv];
    
    
    min_pathdist(j1) = min( abs(P(:,1))+abs(P(:,2))   );  %     COMPUTING THE SHORTEST PATH DISTANCE
    
end

%% THEORETICAL CDF
ll = lambda_l;
lc = lambda_c;
thy_cdf_Tm = @(tm) 1 - exp( -4*lc*tm - 4*ll*tm +2*ll/lc*(1-exp(-2*lc*tm)) );

%% PLOT SIMULATION VS THEORY
resolution = 20;
tm_val = 0:max(min_pathdist)/resolution:max(min_pathdist);
thy_cdf = zeros(size(tm_val)); 
sim_cdf = zeros(size(tm_val));
for k1=1:length(tm_val)
    thy_cdf(k1) = thy_cdf_Tm(tm_val(k1)); % CDF from evaluating the theoretical expression
    sim_cdf(k1) = nnz(min_pathdist<=tm_val(k1))/length(min_pathdist); % CDF from Monte-Carlo simulations
end

figure; hold on
plot(tm_val, thy_cdf, '-r','linewidth',1.5);
plot(tm_val, sim_cdf, 'bo', 'linewidth', 1.5, 'MarkerFaceColor', [0 0 1]);
set(gca, 'FontSize',14)
xlabel('Distance $t_m$ (km)','FontSize',14,'Interpreter','latex');
ylabel('CDF $F_{T_m}(t_m)$','FontSize',14,'Interpreter','latex');
hg=legend('Theorem 1', 'Simulation', 'location','southeast');
set(hg, 'FontSize',14,'Interpreter','latex');
box on

toc