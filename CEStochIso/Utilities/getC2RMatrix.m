function [ C2R, R2C ] = getC2RMatrix( lmax )
%GETC2RMATRIX get the coversion matrix
% from the complex SpH format to the new real SpH format
%
% returns the transformation matrix C2R, and its hermitian transpose R2C
%
% The only non-trivial sub-matrix of C2R is
%               /   (-1)^m      1 \     / +m \     / k=0 \
%  1/sqrt(2) * |                   | * |      | = |       |
%               \ i*(-1)^m     -i /     \ -m /     \ k=1 /
%
% Thus
%    plmk = C2R * plm
% where plmk are the components in the new real SpH format, and
% plm are the components in the complex SpH format.
%
% The new realSpH format has the following ordering:
% 
% index-1  l m k     data(real)
% -----------------------------
%  0       0 0 0       x
%  1       1 0 0       x
%  2       1 1 0       x
%  3       1 1 1       x
%  4       2 0 0       x
%  5       2 1 0       x
%  6       2 1 1       x
%  7       2 2 0       x
%  8       2 2 1       x
%  9       3 0 0       x
% 10       3 1 0       x
% 11       3 1 1       x
% 12       3 2 0       x
% 13       3 2 1       x
% 14       3 3 0       x
% 15       3 3 1       x
% 16       4 0 0       x
% 17       4 1 0       x
% 18       4 1 1       x
% ..       . . .       .
%
% Here l,m are the usual SpH indices, k=0 is the cos component, k=1 is the
% sine component.
%
% This format has the size [(lmax+1)^2,1], i.e. it is a true vector
% i.e. direct vector math is possible, as with the complex SpH format.
%
% Index math:
%       index-1 = l^2 + max(2*m-1,0) + k
% 
%------------------------------------------------------------------
% Contrast this against the old realSpH format (lmcosi -type), that is
% defined by spherelib. It has the following structure: 
%
%   l,m,cosine and sine coefficients of real spherical harmonics,
%   in order of l, then m, i.e.
%
%                       data (4 colums)
% index-1  l m k        l m cos sin
% -------------------------------------
%  0       0 0 0        0 0  x   0
%  1       1 0 0        1 0  x   0
%  2       1 1 [0,1]    1 1  x   x
%  3       2 0 0        2 0  x   0
%  4       2 1 [0,1]    2 1  x   x
%  5       2 2 [0,1]    2 2  x   x
%  6       3 0 0        3 0  x   0
%  7       3 1 [0,1]    3 1  x   x
%  8       3 2 [0,1]    3 2  x   x
%  9       3 3 [0,1]    3 3  x   x
% 10       4 0 0        4 0  x   0
% 11       4 1 [0,1]    4 1  x   x
% ..       . . .        . .  .   .
%
% Index math:
%       index-1 = l*(l+1) + m
%
%------------------------------------------------------------------
%
% Written by Stefan Ballmer
% sballmer@ligo.caltech.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N=(lmax+1)^2;
% initialize the index vectors for the new realSpH format
    indM1=[0:N-1]';
    lrvec=floor(sqrt(indM1));
    mrvec=floor((indM1+1-lrvec.^2)/2);
    krvec=indM1 - lrvec.^2 - max(2*mrvec-1,0);
    
    lrM= lrvec*ones(1,N);
    mrM= mrvec*ones(1,N);
    krM= krvec*ones(1,N);

% initialize the index vectors for the complex SpH format
    [lvec,mvec]=getLMvec(lmax);
    lM = ones(N,1)*lvec';
    mM = ones(N,1)*mvec';
    
% build the C2R matrix
    M=zeros(N,N);
    
    % if m=0: trivial transform
    ind= and( lrM==lM , and( mrM==0 , mM==0 ) );
    M(ind)=1;
    
    % m positive, k=0: (-1)^m / sqrt(2)
    ind= and( lrM==lM , and( mrM~=0 , and( mrM==mM , krM==0 )));
    M(ind)=    (-1).^mrM(ind) / sqrt(2);

    % m positive, k=1: 1i*(-1)^m / sqrt(2)
    ind= and( lrM==lM , and( mrM~=0 , and( mrM==mM , krM==1 )));
    M(ind)=1i* (-1).^mrM(ind) / sqrt(2);

    % m negative, k=0: 1 / sqrt(2)
    ind= and( lrM==lM , and( mrM~=0 , and( mrM==-mM , krM==0 )));
    M(ind)=   1 / sqrt(2);

    % m negative, k=1: -1i / sqrt(2)
    ind= and( lrM==lM , and( mrM~=0 , and( mrM==-mM , krM==1 )));
    M(ind)= -1i / sqrt(2);

    % only identical l and abs(m) are non-zero
    ind= or( lrM~=lM , mrM~=abs(mM) );
    M(ind)=0;
    
    C2R = M;
    if nargout>1
        R2C = M'; % M is unitary
    end
end

