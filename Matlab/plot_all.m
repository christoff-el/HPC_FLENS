function plot_all(procCount)

figure;
hold on

for i=0:procCount-1
    
    plot_solution(i);
    
end

view(64.5,22);
