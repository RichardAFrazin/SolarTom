function y = Aregmult(x,mode,A,Rs)

    if isa(A,'numeric')
        explicitA = true;
    elseif isa(A,'function_handle')
        explicitA = false;
    else
        error('Aregmult:','%s','A must be numeric or a function handle');
    end
    
    [m,n] = size(A);
    [s,n1] = size(Rs);
    
    if (n1 ~= n)
      error(['Aregmult: Rs and A must have the same number of ' ...
             'columns!']);
    end
    
    switch mode
        case 0     % return size
            y=[m+s,n];
        case 1      % y = A*x
            if explicitA
                y = [ A*x ; Rs*x ];
            else
                y = [ A(x,1) ; Rs*x ];
            end
        case 2      % y = A'*x
            if explicitA
                y = A'*x(1:m) + Rs'*x(m+1:m+s) ;
            else
                y = A(x(1:m),2) + Rs'*x(m+1:m+s) ;
            end        
    end

    return
    end