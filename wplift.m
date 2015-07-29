function [C,y2,dec]=wplift(xa,U,P,level)
%@Wavelet packet decomposition tree Ruben Dario Pinzon Morales 2010
%E : catracteristics including normalized energy
%y : decomposition without filter mirror correction
%y2: decomposition with filtrer mirror correction
%dec:decomposition indexes 

y={xa};
index=[2 3];
base=1;
for j=1:level
    for node=1:length(base)
        [y{index(1+2*(node-1))},y{index(2*(node))}]=lift(y{base(node)},U,P,1);
    end
    base=index;
    index=index(end)+1:index(end)+2^(j+1);
end
%Filter mirror correction
G=1;
index=1;
dec=zeros(1,level);
dec(1)=length(xa);
y2=y;
for j=1:level
    G = grays(j)+1+max(G);
    for sub=1:length(G)
        y2{index+sub}=y{G(sub)};
        lel(index+sub)=j;
    end
    index=max(G);
    dec(j+1)=length(xa)/2^(j);
end

%% ============= Normalized linear Energy computing for best basis selection
E=zeros(1,length(y));
for nodes=1:length(y)
    E(nodes)=log(sum(y{nodes}.^2)/(length(y{nodes})/(2^lel(nodes))));
    %E(nodes)=sum(y{nodes}.^2)/length(y{nodes});
    %car(nodes,:)=measures(y{nodes});
end
E=E./E(1);
%C=[E' car];
C=E;





