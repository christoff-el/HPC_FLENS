function show_mesh_sol(elements, coordinates,  solution , varargin)
  art=['r*'; 'bx';'mo'];

 % triplot(elements',coordinates(1,:)',coordinates(2,:)', 'Color', 'g');
% figure;
  trisurf(elements,coordinates(:,1),coordinates(:,2),solution,'facecolor','interp','edgecolor','none')
%   figure
%   [isolinex,isoliney] = isolines(elements,coordinates,solution(elements),100);
%   plot(isolinex',isoliney','k-','linewidth',0.8);
%   axis tight
  
  
%   for j=1:nargin-3;
%     Boundary=varargin{j};
%     for i=1:size(Boundary,2)
%        line(coordinates(1,Boundary(:,i)), coordinates(2,Boundary(:,i)),'Color',art(j,1));
%     end
%   end
%
%   center=zeros(2,size(elements,2));
%   for i=1:size(elements,2)
%     center(:,i)=sum(coordinates(:,elements(:,i)),2)/3;
%   end
%
%   elementnum=1:1:size(elements,2);
%   text(center(1,:),center(2,:),num2str(elementnum'), 'color', 'red', ...
%       'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'red' );
%
%   koordnum=1:1:size(coordinates,2);
%   text(coordinates(1,:),coordinates(2,:), num2str(koordnum'),'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

