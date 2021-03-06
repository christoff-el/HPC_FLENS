serial_cg_counts = [33025 131585 525313 2099201];
serial_assembly = [0.817531 3.41866 14.7845 65.8931];
serial_cg = [1.49912 12.4984 101.033 802.526];
serial_gs_counts = 33025;
serial_gs = 951.417;

serial_goto_assembly = [0.827823 3.54683 15.4809 67.2637];
serial_goto_cg = [1.02973 9.27832 78.1647 620.832];
serial_goto_gs = 958.206;

par_cg_counts = [33284 132100 526340 2101252 8396804];
par_assembly = [0.195 0.824 3.57 15.27 65];
par_cg = [0.8422 6.5692 54.327 459.899 3641.48];
par_gs_counts = [8452 33284];
par_gs = [12.11 184.962];

par_goto_assembly = [0.204 0.845 3.68 16.01 68.1];
par_goto_cg = [0.551 4.468 40.657 329.28 2589.78];

lw=1.3;
figure; 
loglog(serial_cg_counts, serial_assembly,'-rx','LineWidth',2*lw); hold on;
loglog(serial_cg_counts, serial_goto_assembly,'-gx','LineWidth',lw);
loglog(par_cg_counts, par_assembly, '-bx','LineWidth',2*lw);
loglog(par_cg_counts, par_goto_assembly, '-mx','LineWidth',lw);
xlabel('Mesh Element Count');
ylabel('Computation time');
title('Serial vs. Parallel for System Assembly')
legend('Serial','Serial with GotoBLAS','Parallel','Parallel with GotoBLAS','Location','NorthWest')

matlab2tikz('assemblyPlot.tikz', 'height', '\figureheight', 'width', '\figurewidth','strict',true);
close all;

figure; 
loglog(serial_cg_counts, serial_cg,'-rx','LineWidth',2*lw); hold on;
loglog(serial_cg_counts, serial_goto_cg,'-gx','LineWidth',lw);
loglog(par_cg_counts, par_cg, '-bx','LineWidth',2*lw);
loglog(par_cg_counts, par_goto_cg, '-mx','LineWidth',lw);
xlabel('Mesh Element Count');
ylabel('Computation time');
title('Serial vs. Parallel for System Solving')
legend('Serial CG','Serial CG with GotoBLAS','Parallel CG','Parallel CG with GotoBLAS','Location','NorthWest')
matlab2tikz('solvePlot.tikz', 'height', '\figureheight', 'width', '\figurewidth','strict',true);
close all;

