function op = linop_TVthetaphi4D( sz, action )

%LINOP_TVthetaphi   TFOCS 2D Total-Variation (TV) linear operator
%             on theta,phi for a spherical object 
%
%   inspired from LINOP_TV TFOCS v1.1 
%    D.VIBERT 02/2012
%
%    OP = LINOP_TVthetaphi( SZ ) returns a handle to a TFOCS linear operator that
%      implements the total variation linear operator on an 
%      Mrad x Ntheta x Ophi grid in all the theta,phi plane only
%      that is, to be applied to 3D spherical objects of size [Mrad,Ntheta,Ophi].
%      By default, it expects to operate on M*N*O x 1 vectors
%      but if SZ = {M,N,O}, then expects to operate on M x N x O matrices
%
%    TV = LINOP_TVthetaphi( X ) returns ||X||_TV  if X is bigger than 2 x 2
%
%    [...] = LINOP_TV(SZ, ACTION )
%       if ACTION is 'handle', returns a TFOCS function handle (default)
%       if ACTION is 'cvx', returns a function handle suitable for CVX
%       if ACTION is 'matrix', returns the explicit TV matrix
%           (real part corresponds to horizontal differences,
%            imaginary part correspond to vertical differences)
%       if ACTION is 'norm', returns an estimate of the norm
%
%   TODO: check circular case for non-square domains
%

error(nargchk(1,2,nargin));
if nargin < 2 || isempty(action), action = 'handle'; end

CALCULATE_TV = false;
if numel(sz) > 5 
    CALCULATE_TV = true;
    X   = sz;
    sz  = size(X);
end

if iscell(sz)
    n1 = sz{1}; % nrad
    n2 = sz{2}; % ntheta
    n3 = sz{3}; % nphi
    n4 = sz{4}; % ntime
else
    n1 = sz(1); % nrad
    n2 = sz(2); % ntheta
    n3 = sz(3); % nphi
    n4 = sz(4); % ntime
end

if strcmpi(action,'norm')
    n1 = 1; %to reduce norm computation, since independant from n1
    n4 = 1;
end

% Setup the Total-Variation operators
mat = @(x) reshape(x,n1,n2,n3,n4);

if strcmpi(action,'matrix') || strcmpi(action,'cvx')
    e = ones(max(n2,n3),1);
    e2 = e;
    J = spdiags([-e2,e], 0:1,n3,n3);
    J(end,1) = 1; % circular in phi
    I = speye(n2);
    Dh = kron(J,I); 
    Dh = kron(Dh,speye(n1)); % phi differences, sparse matrix
    Dh = kron(speye(n4),Dh);
   
    
    e2 = e;
    e2(n2:end) = 0;
    J = spdiags([-e2,e], 0:1,n2,n2);
    I = speye(n3);
    Dv = kron(I,J);  
    Dv = kron(Dv,speye(n1)); % theta differences, sparse matrix
    Dv = kron(speye(n4),Dv);
    
    if strcmpi(action,'matrix')
        op = Dh + 1i*Dv;
    else
        % "norms" is a CVX function
        op = @(X) sum( norms( [Dh*X(:), Dv*X(:)]' ) );
    end
    return;
end

%Dphi , circular
Dh     = @(X) vec( cat(3, diff(X,1,3), X(:,:,1,:) - X(:,:,end,:)) );
diff_h = @(X) cat(3, X(:,:,end,:), X(:,:,1:end-1,:)) - X;
%Dtheta not circular
Dv     = @(X) vec( cat(2, diff(X,1,2), zeros(n1,1,n3,n4)) );
diff_v = @(X) cat(2, zeros(n1,1,n3,n4), X(:,1:end-1,:,:))...
        - cat(2, X(:,1:end-1,:,:), zeros(n1,1,n3,n4));

if iscell(sz)
    Dh_transpose = @(X)      diff_h(mat(X))  ;
    Dv_transpose = @(X)      diff_v(mat(X))  ;
else
    Dh_transpose = @(X) vec( diff_h(mat(X)) );
    Dv_transpose = @(X) vec( diff_v(mat(X)) );
end

TV  = @(x) ( Dh(mat(x)) + 1i*Dv(mat(x)) );     % real to complex
TVt = @(z) ( Dh_transpose(real(z)) + Dv_transpose(imag(z)) );

if CALCULATE_TV
    op = norm( TV(X), 1 );
    return;
end

if strcmpi(action,'norm')
    % to compute max eigenvalue, I use a vector
    % that is very likely to be the max eigenvector:
    %  matrix with every entry alternating -1 and 1
    even = @(n) ~( n - 2*round(n/2) );  % returns 1 if even, 0 if odd
    Y = zeros( n2 + even(n2), n3 + even(n3));
    nn = numel(Y);
    Y(:) = (-1).^(1:nn);
    Y = Y(1:n2,1:n3);
    Y = reshape(Y,1,n2,n3,1);
    op = norm( TV(Y) )/norm(Y(:));
    
    % Nearly equivalent to:
    % norm(full( [real(tv); imag(tv)] ) ) 
    % where tv is the matrix form
else
    if iscell(sz)
        szW = { [n1,n2,n3,n4], [n1*n2*n3*n4,1] };
    else
        szW = n1 * n2 * n3 * n4;
        szW = [szW,szW];
    end
    op = @(x,mode)linop_tv_r2c(szW,TV,TVt,x,mode);
end

function y = linop_tv_r2c( sz, TV, TVt, x, mode )
switch mode,
    case 0, y = sz;
    case 1, y = TV( realcheck( x ) );
    case 2, y = realcheck( TVt( x ) );
end

function y = realcheck( y )
if ~isreal( y ), 
    error( 'Unexpected complex value in linear operation.' );
end

