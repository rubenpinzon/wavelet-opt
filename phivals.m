function [x,phi,psi] = phivals(h,i)
% function [x,phi,psi] = phivals(h,i)
% Generate a scaling function and its associated wavelet
% using the given filter coefficients
%               Kevin Amaratunga
%               5 March, 1993
%
%     h = filter coefficients (sum(h)=2) 
%     i = discretization parameter.  The number of points per integer
%         step is 2^i.  Thus, setting i = 0 gives the scaling function
%         and wavelet values at integer points.
%
if i < 0
  error('phivals: i must be non-negative')
end
[m,n] = size(h);
g = (-1).^(0:m-1)' .* h(m:-1:1);

% The Haar filter produces a singular matrix, but since we know the solution
% already we treat this as a special case.
if m == 2 & h == [1;1]
  phi = [ones(2^i,1);0];
  if i > 0
    psi = [ones(2^(i-1),1);-ones(2^(i-1),1);0];
  elseif i == 0
    psi = [1;0];
  end
else
  ch = [h; zeros(m,1)];
  rh = [h(1), zeros(1,m-1)];
  tmp = toeplitz(ch,rh);
  M = zeros(m,m);
  M(:) = tmp(1:2:2*m*m-1);
  M = M - eye(m);
  M(m,:) = ones(1,m);
  tmp = [zeros(m-1,1); 1];
  phi = M \ tmp;		% Integer values of phi

  if i > 0
    for k = 0:i-1
      p = 2^(k+1) * (m-1) + 1;	% No of rows in toeplitz matrix
      q = 2^k * (m-1) + 1;	% No of columns toeplitz matrix
      if (k == 0)
        ch0 = [h; zeros(p-1-m,1)];
        ch = [ch0; 0];
        cg0 = [g; zeros(p-1-m,1)];
      else
        ch = zeros(p-1,1);
        ch(:) = [1; zeros(2^k-1,1)] * ch0';
        ch = [ch; 0];
      end
      rh = [ch(1), zeros(1,q-1)];
      Th = toeplitz(ch,rh);
      if k == i-1
        cg = [1; zeros(2^k-1,1)] * cg0';
	cg = cg(:);
        cg = [cg; 0];
        rg = [cg(1), zeros(1,q-1)];
        Tg = toeplitz(cg,rg);
        psi = Tg * phi;
      end
      phi = Th * phi;
    end
  elseif i == 0
    cg0 = [g; zeros(m-2,1)];
    cg = [cg0; 0];
    rg = [cg(1), zeros(1,m-1)];
    Tg = toeplitz(cg,rg);
    psi = Tg * phi;
    psi = psi(1:2:2*m-1);
  end
end

[a,b] = size(phi);
x = (0:a-1)' / 2^i;
