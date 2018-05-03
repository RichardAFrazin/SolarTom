function op = linop_dt4D( sz, action )

%LINOP_dt4D   TFOCS  time derivative linear operator
%                      on a 4D (rad,theta,phi,time) spherical object
%
%   inspired from LINOP_TV TFOCS v1.1 
%    D.VIBERT 05/2012
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
if numel(sz) > 4 
    CALCULATE_TV = true;
    X   = sz;
    sz  = size(X);
end

if iscell(sz)
    n1 = sz{1}; % nspace
    n2 = sz{2}; % ntime
else
    n1 = sz(1); % nspace
    n2 = sz(2); % ntime
end

if strcmpi(action,'norm')
    n1 = 1; %to reduce norm computation, since independant from n1
end

% Setup the Total-Variation operators
mat = @(x) reshape(x,n1,n2);

if strcmpi(action,'matrix') || strcmpi(action,'cvx')
    e = ones(n2,1);
    e2 = e;
    J = spdiags([-e2,e], 0:1,n2,n2);
    I = speye(n1);
    Dt = kron(J,I); % time differences, sparse matrix
    
    if strcmpi(action,'matrix')
        op = Dt;
    else
        % "norms" is a CVX function
        op = @(X) sum( abs( Dt*X(:) ) ); 
    end
    return;
end

%Dt
Dt     = @(X) vec( cat(2, diff(X,1,2), zeros(n1,1)) );
diff_t = @(X) cat(2, zeros(n1,1), X(:,1:end-1)) - cat(2, X(:,1:end-1), zeros(n1,1));

DT4D  = @(x) ( Dt(mat(x)) );
if iscell(sz)
    DT4Dt = @(X)      diff_t(mat(X))  ;
else
    DT4Dt = @(X) vec( diff_t(mat(X)) );
end

if CALCULATE_TV
    op = norm( DT4D(X), 1 );
    return;
end

if strcmpi(action,'norm')
    % to compute max eigenvalue, I use a vector
    % that is very likely to be the max eigenvector:
    %  matrix with every entry alternating -1 and 1
    Y = (-1).^(1:n2);
    Y = reshape(Y,1,n2);
    op = norm( DT4D(Y) )/norm(Y(:));
else
    if iscell(sz)
        szW = { [n1,n2], [n1*n2,1] };
    else
        szW = n1 * n2;
        szW = [szW,szW];
    end
    op = @(x,mode)linop_dt4D_mode(szW,DT4D,DT4Dt,x,mode);
end

function y = linop_dt4D_mode( sz, DT4D, DT4Dt, x, mode )
switch mode,
    case 0, y = sz;
    case 1, y = DT4D( x );
    case 2, y = DT4Dt( x ) ;
end


