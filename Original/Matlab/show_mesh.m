function show_mesh(elements, coordinates,proc,flag,show_text,varargin)
  art=['r*'; 'bx';'mo'];
  %*** define colormap
  if flag==1
    trisurf(elements,coordinates(:,1),coordinates(:,2),proc*ones(size(coordinates,1),1)');
    view(2);
  else
    triplot(elements,coordinates(:,1),coordinates(:,2), 'Color', 'g');
    for j=1:nargin-5;
      Boundary=varargin{j};
      for i=1:size(Boundary,1)
         line(coordinates(Boundary(i,:),1), coordinates(Boundary(i,:),2),'Color',art(j));
      end
    end
    center=zeros(size(elements,1),2);
    for i=1:size(elements,1)
      center(i,:)=sum(coordinates(elements(i,:),:),1)/3;
    end
  
    if show_text 
      elementnum=1:size(elements,1);
      text(center(:,1),center(:,2),num2str(elementnum'), 'color', 'red', ...
          'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'red' );

      koordnum=1:size(coordinates,1);
      text(coordinates(:,1),coordinates(:,2), num2str(koordnum'),'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
  end
  
  axis equal;
