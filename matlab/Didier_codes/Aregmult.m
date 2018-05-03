function y = Aregmult(x,mode,A,Rs,Rt)

    if isa(A,'numeric')
        explicitA = true;
    elseif isa(A,'function_handle')
        explicitA = false;
    else
        error('Aregmult:','%s','A must be numeric or a function handle');
    end
    
    [m,n] = size(A);
    [s,n] = size(Rs);
    [t,n] = size(Rt);
    
    switch mode
        case 0     % return size
            y=[m+s+t,n];
        case 1      % y = A'*x
            if explicitA
                y = [ A*x ; Rs*x ; Rt*x ];
            else
                y = [ A(x,1) ; Rs*x ; Rt*x ];
            end
        case 2      % y = A*x
            if explicitA
                y = A'*x(1:m) + Rs'*x(m+1:m+s) + Rt'*x(m+s+1:m+s+t);
            else
                y = A(x(1:m),2) + Rs'*x(m+1:m+s) + Rt'*x(m+s+1:m+s+t);
            end        
    end
    
%     if strcmp(transp_flag,'notransp')      % y = A'*x
%         if explicitA
%             y = [ A*x ; lambda*(Rs*x) ; mu*(Rt*x) ];
%         else
%             y = [ A(x,'notransp') ; lambda*(Rs*x) ; mu*(Rt*x) ];
%         end
%     elseif strcmp(transp_flag,'transp') % y = A*x
%         if explicitA
%             y = A'*x(1:m)+lambda*(Rs'*x(m+1:m+s))+mu*(Rt'*x(m+s+1:m+s+t));
%         else
%             y = A(x(1:m),'transp')+lambda*(Rs'*x(m+1:m+s))+mu*(Rt'*x(m+s+1:m+s+t));
%         end
%     end    
%     
    return
    end