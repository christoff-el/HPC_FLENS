function show_edgeno(coordinates,varargin)

  art=['r';'b'; 'm'];
  for j=1:nargin-1;
    edge2nodes=varargin{j};
    edgeno=1:1:size(edge2nodes,2);
    XBmid=(coordinates(1,edge2nodes(1,:))+coordinates(1,edge2nodes(2,:)))/2;
    YBmid=(coordinates(2,edge2nodes(1,:))+coordinates(2,edge2nodes(2,:)))/2;
    text(XBmid,YBmid,num2str(edgeno'),...
       'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','BackgroundColor',art(j,1));
  end