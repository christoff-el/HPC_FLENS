serial_cg_counts = [33025 131585 525313 2099201];
serial_assembly = [0.157823 0.710118 2.9903 14.174];
serial_cg = [0.520034 5.3675 44.6949 363.397];

serial_goto_assembly = [0.827823 3.54683 15.4809 67.2637];
serial_goto_cg = [1.02973 9.27832 78.1647 620.832];

par_cg_counts = [33284 132100 526340 2101252 8396804];
par_assembly = [0.195 0.824 3.57 15.27 65];
par_cg = [0.8422 6.5692 54.327 459.899 3641.48];

par_goto_assembly = [0.204 0.845 3.68 16.01 68.1];
par_goto_cg = [0.551 4.468 40.657 329.28 2589.78];

orig_serial_cg_counts = [33025 131585 525313 2099201];
orig_serial_assembly = [0.080795 0.328981 1.30698 ];
orig_serial_cg = [0.514473 5.24624 44.4589 ];

orig_par_cg_counts = [33284 132100 526340 2101252 8396804];
orig_par_assembly = [];
orig_par_cg = [];

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
legend('Serial CG','Serial CG with GotoBLAS','Parallel CG','Parallel G with GotoBLAS','Location','NorthWest')
matlab2tikz('solvePlot.tikz', 'height', '\figureheight', 'width', '\figurewidth','strict',true);
close all;

