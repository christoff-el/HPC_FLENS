function [] = plot_triangulation(proc,flag, show_text)


if proc == 0
  coord_str = sprintf('../examples/output/coordinates.dat');
  elements_str = sprintf('../examples/output/elements.dat');
  dirichlet_str = sprintf('../examples/output/dirichlet.dat');
  neumann_str = sprintf('../examples/output/neumann.dat');
  % edge2nodes_str = sprintf('../examples/output/edge2nodes.dat');
else
  coord_str = sprintf('../examples/output/coordinates%d.dat', proc);
  elements_str = sprintf('../examples/output/elements%d.dat', proc);
  dirichlet_str = sprintf('../examples/output/dirichlet%d.dat', proc);
  neumann_str = sprintf('../examples/output/neumann%d.dat', proc);
  % edge2nodes_str = sprintf('../examples/output/edge2nodes%d.dat', proc);
end

%*** load mesh files
coordinates = load(coord_str);
elements = load(elements_str);
dirichlet = load(dirichlet_str);
neumann = load(neumann_str);
% edge2nodes = load(edge2nodes_str);

% Show Mesh with coloured boundary edges
show_mesh(elements, coordinates, proc,flag,show_text, dirichlet,neumann);
%   -> Inner edges green
%   -> neumann edges blue
%   -> dirichlet edges red

% Plot Edge numbers
if flag==0
%   show_edgeno(coordinates',edge2nodes');
end
