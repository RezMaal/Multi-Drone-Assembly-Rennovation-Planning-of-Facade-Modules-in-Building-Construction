function T=ResTable(Route,Rb)
l=size(Route,1);
Tc=sum(min(Rb.C,[],2))/Rb.m;
T=zeros(l,4);
for i=1:l
    T(i,:)=[max(Route(i,:)),min(Route(i,:)),mean(Route(i,:)),mean(abs(Route(i,:)-Tc))];
end
end