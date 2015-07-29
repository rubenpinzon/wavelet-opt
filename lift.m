function [CA,CD,index]=lift(y,U,P,level)
%LIFT    Analisys based on lifting scheme 
%
%
%
%Ruben Dario Pinzon Morales @2009
if size(y,2)==1
    y=y';
end
lo=length(y);                       %Data length
xe=y(1:2:lo);                       %Even samples
xo=y(2:2:lo);                       %Odd samples
N=length(U);                        %U Filter order
N_p=length(P);                      %P Filter order
n=length(xe);                       %signal's length

Xe=xe;
Xo=xo;
for p=1:floor(N_p/2)                  %signal padding
    Xe=[xe(p) Xe xe(end-p)];
    Xo=[xo(p) Xo xo(end-p)];
end
CD=zeros(1,floor(lo/2));
for i=floor(N_p/2)+1:n+floor(N_p/2)
    CD(i-floor(N_p/2))=Xo(i)-P*Xe(i-floor(N_p/2):i+floor(N_p/2)-1)';    
end

cd=CD;
for p=1:floor(N/2)                  %signal padding
    cd=[CD(p) cd CD(end-p)];
end
index=n;
CA=zeros(1,floor(lo/2));
for i=floor(N/2)+1:n+floor(N/2)
    CA(i-floor(N/2))=Xe(i)+U*cd(i-floor(N/2)+1:i+floor(N/2))';    
end

x=CA;
if level>1                          % Recursion if level > 1.
   [ca,cd,i]=lift(x,U,P,level-1);
   index=[index i];
   CA=[CA ca];
   CD=[CD cd];
   %x=ca;
   clear ca cd
end
