%%% Haswell Architecture
beta_l1 = 12*8;
beta_l2 = 8*8;
beta_l3 = 4*8;
beta_ram = 2*8;
pi_scalar_fma = 2*2;  % 2 execution units w/ FMA
pi_simd = 2*8;   % 2 execution units w/ FMA using single precision float
pi_simd_fma = 2*2*8;   % 2 execution units w/ FMA using single precision float
pi_simd_reduced = pi_simd*3/4*35/40;
pi_simd_fma_reduced = pi_simd_fma*3/4*35/40;
% TODO add latency bounded peak performance of 11.62 ish

%%% Plot Bounds (as powers of 2)
I_bounds = [-7 4];
P_bounds = [-2 6];

%%% Colors
k = [0.34 0.34 0.34];

%%% Testing parameters and results
window_size = 35;

W_exp = 0;
W_baseline =                       48*window_size^2 + W_exp*window_size^2;
W_rgba     = 19 +  4*window_size + 38*window_size^2;
W_fastpre  =  3 +  4*window_size + 37*window_size^2;
W_simd_v3  = 43 + 52*window_size + 49*window_size^2;
W_simd_v7  = 43 + 52*window_size + 49*window_size^2;
W_simd_v10 = 62 + 49*window_size + 51*window_size^2;
W_simd_v16 = 62 + 49*window_size + 51*window_size^2;
W_simd_v17 = 62 + 49*window_size + 51*window_size^2;

Q_baseline =                      29*window_size^2;
Q_rgba     = 12 + 0*window_size + 25*window_size^2;
Q_fastpre  = 12 + 0*window_size + 25*window_size^2;
Q_simd_v3  = 12 + 0*window_size + 28*window_size^2;
Q_simd_v7  = 12 + 0*window_size + 28*window_size^2;
Q_simd_v10 = 12 + 0*window_size + 68*window_size^2;
Q_simd_v16 = 12 + 0*window_size + 68*window_size^2;
Q_simd_v17 = 12 + 0*window_size + 68*window_size^2;

I_baseline = W_baseline/Q_baseline;
I_rgba     = W_rgba/Q_rgba;
I_fastpre  = W_fastpre/Q_fastpre;
I_simd_v3  = W_simd_v3/Q_simd_v3;
I_simd_v7  = W_simd_v7/Q_simd_v7;
I_simd_v10 = W_simd_v10/Q_simd_v10;
I_simd_v16 = W_simd_v16/Q_simd_v16;
I_simd_v17 = W_simd_v17/Q_simd_v17;

P_baseline = W_baseline/110129;
P_rgba     = W_rgba/48893;
P_fastpre  = W_fastpre/40031;
P_simd_v3  = W_simd_v3/16370;
P_simd_v7  = W_simd_v7/16540;
P_simd_v10 = W_simd_v10/12980;
P_simd_v16 = W_simd_v16/12610;
P_simd_v17 = W_simd_v17/11633;

figure();
hold on;

% % Show the performance of all 4 functions (first for easy legend)
scatter(log2(I_baseline), log2(P_baseline), 'ko');
% scatter(log2(I_rgba),     log2(P_rgba),     'k^');
scatter(log2(I_fastpre),  log2(P_fastpre),  'ks');
scatter(log2(I_simd_v3),  log2(P_simd_v3),  'kx');
% scatter(log2(I_simd_v7),  log2(P_simd_v7),  'k+');
scatter(log2(I_simd_v10), log2(P_simd_v10), 'kp');
% scatter(log2(I_simd_v16), log2(P_simd_v16), 'k<');
scatter(log2(I_simd_v17), log2(P_simd_v17), 'k<');

text(log2(I_baseline)+0.3, log2(P_baseline)     , 'baseline'          , 'HorizontalAlign','left'); 
% text(log2(I_rgba)    -0.3,     log2(P_rgba)     , 'rgba'              , 'HorizontalAlign','right'); 
text(log2(I_fastpre) +0.3,  log2(P_fastpre)+.35 , 'fastexp+'          , 'HorizontalAlign','left');
text(log2(I_fastpre) +0.3,  log2(P_fastpre)+.05 , 'precompute'        , 'HorizontalAlign','left');
text(log2(I_simd_v3) +1.5,  log2(P_simd_v3)-.25 , 'simd v3'           , 'HorizontalAlign','right'); 
% text(log2(I_simd_v7) +0.3,  log2(P_simd_v7)     , 'simd v7'           , 'HorizontalAlign','left'); 
text(log2(I_simd_v10)-0.3, log2(P_simd_v10)-.15 , 'simd v10'          , 'HorizontalAlign','right'); 
% text(log2(I_simd_v16)-0.3, log2(P_simd_v16)+.15 , 'simd v16'          , 'HorizontalAlign','right'); 
text(log2(I_simd_v17)+0.3, log2(P_simd_v17)+.15 , 'simd v17'          , 'HorizontalAlign','left'); 

% Plot beta lines for different cache levels
plot(I_bounds, log2(beta_l1*2.^(I_bounds)), 'k');
plot(I_bounds, log2(beta_l2*2.^(I_bounds)), 'k');
plot(I_bounds, log2(beta_l3*2.^(I_bounds)), 'k');
plot(I_bounds, log2(beta_ram*2.^(I_bounds)), 'k:');

beta_ang = 45;
text(I_bounds(1)+.25,  0.20, "\beta_{L1}=96", 'rotation',beta_ang);
text(I_bounds(1)+.50, -0.9 , "\beta_{L2}=64 bytes/cycle", 'rotation',beta_ang);
text(I_bounds(1)+.75, -1.75, "\beta_{L3}=32 bytes/cycle", 'rotation',beta_ang);
text(I_bounds(1)+1.75, -1.75, "\beta_{ram}=16 bytes/cycle", 'rotation',beta_ang);

% Plot peak performance single core
plot(I_bounds, log2([pi_scalar_fma,pi_scalar_fma]), 'k');
plot(I_bounds, log2([pi_simd,pi_simd]), 'k');
plot(I_bounds, log2([pi_simd_fma,pi_simd_fma]), 'k');
plot(I_bounds, log2([pi_simd_reduced,pi_simd_reduced]), 'k:');
plot(I_bounds, log2([pi_simd_fma_reduced,pi_simd_fma_reduced]), 'k:');

text(I_bounds(1)+0.25, log2(pi_scalar_fma)+.6, "Scalar FMA:");
text(I_bounds(1)+0.25, log2(pi_scalar_fma)+.3, "\pi=4");
text(I_bounds(1)+0.25, log2(pi_simd)+.25, "SIMD: \pi=16 flops/cycle");
text(I_bounds(1)+0.25, log2(pi_simd_fma)+.25, "SIMD FMA: \pi=32 flops/cycle");
text(I_bounds(2)-0.25, log2(pi_simd_reduced)+.25, "\pi=10.5 flops/cycle", 'HorizontalAlign','right');
text(I_bounds(2)-0.25, log2(pi_simd_fma_reduced)+.25, "\pi=21 flops/cycle", 'HorizontalAlign','right');

% Setup the rest of the plot
axis equal;
axis([min(I_bounds) max(I_bounds) min(P_bounds) max(P_bounds)]);
grid on;
set(gca,'GridColor',[1,1,1]);
set(gca,'GridAlpha',1);
title("Haswell CPU (Intel(R) Xeon(R) CPU E5-2680 v3)");
xlabel("Operational Intensity (flops/byte)");
ylabel("Performance (flops/cycle)");
% legend("P_{mcost}", "P_{spatial}", "P_{view}", "P_{plane}", "P_{total}", 'Location', 'South', 'NumColumns',5);
xticklabels(arrayfun(@(x) sprintf('2^{%d}',x), xticks, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('2^{%d}',y), yticks, 'UniformOutput', false));
set(gca,'Color',[.92,.92,.92])
set(gcf,'color','w');
box off
