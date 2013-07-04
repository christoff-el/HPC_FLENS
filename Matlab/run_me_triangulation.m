clear all; close all; clc;

% find out how many output files there are
proc = check_num_files('../examples/output/coordinates')+1;

show_text = 1;

if proc>1
    % plot mesh for each subdomain
    for i=1:proc
        figure(i)
        plot_triangulation(i-1,0, show_text);
    end
    % plot global mesh
    figure
    hold on
    for i=1:proc
      plot_triangulation(i-1,1, show_text);
    end
    tmp1 = [(1:proc)./(proc*2)+0.4]';
    colormap(cool);
    brighten(0.5);
    
    hold off
else
    figure
    plot_triangulation(0,0,show_text);
end