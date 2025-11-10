function [eig_bounds, ind_range] = veig(A, B, ind)
%%  VEIG: Compute eigenvalue bounds for problem Ax=lambda Bx.
%  Input: 
%     A, B: symmetric real or interval matrix with size n;  B: positive definite. 
%     ind: the eigenvalue index or continuous indices of eigenvalues wanted;.
%  Output:
%     eig_bounds: eigenvalue bounds in interval type;
%     ind_range : indices corresponding to returned eigenvalues.
%
%  Notice:
%  1) The returned result may have more eigenvalues than the ones specified by input parameter "ind".
%  2) This function is not expected to solve large scaled problem even the problem is a sparse one.
%  Reason: The block ldl faction method has been used to verify the bound of eigenvalues.
%  Although the block LDL algorithm can solve sparse problem efficiently, 
%  the calling of standard EIG function to compute ALL approximate eigenvalues
%  may take long time.
%
%  Reference: 
%  1) Yamamoto N., "A simple method for error bounds of eigenvalues of 
%  symmetric matrices", Linear Algebra and its Application, 234(2001), 227-234
%
%  Xuefeng Liu, xfliu.math@gmail.com
%
%  2011/03/11 First version 
%  2011/09/30 Second version. Document modified.
%  2011/11/22 Third version. Modified to give exact bound if possible. Especially the part for repeated eigenvalues.
%  2025/11/7  Revision: Change the proces to estimate the ldl residual error from very small initial value.

%%
eig_bounds = [];
ind_range = [];
%Check the symmetricity of matrix 
if ~isequal(A,A')
  error('Matrix A is not symmetric.')
end

if ~isequal(B,B')
  error('Matrix B is not symmetric.')
end

if size(A,1)>1000
    error("The function veig is designed for dense matrix with small scale. Make sure you select the proper eigenvalue solver.")
end

if ~exist("ind","var")
    min_ind=1; max_ind=size(A,1);
else
    if ~isnumeric(ind)
      error('The 3rd input parameter should be numeric.')
    end
    min_ind=min(ind); max_ind=max(ind);
end

if issparse(A) || issparse(B)
    eig_list = eig(full(mid(A)), full(mid(B)));
else
    eig_list = eig(mid(A), mid(B));
end

eig_bounds = intval([]);

min_eig_B=min(eig(mid(B)));

index=min_ind;

while index <= max_ind
    
    eig_test = eig_list(index);
    
    % Find the possible enlarged range of eigenvalue indices by seeking 
    % rough lower bound and upper bound for the specified index-th eigenvalue.

    lambda_test = eig_test;
    find_lambda = false; k = 0;
    if index == 1
       neg_num_low = 0;
       lambda_test_low = eig_test;
    else
        while ~find_lambda
          %Smallest relative perturbation is about 10^(-10).
          if(k>10)
                error('Too many loops for finding rough lower bound.')
          end
          E = mid(A) - lambda_test*mid(B);
            
          [L,D,P] = ldl(E);
                
          [negative_num, positive_num, zero_num, F] = GetInertia(D);
          if negative_num >= index || F ~= 0
            %Reduce the value of eig_est by relative magnitude.
            lambda_test = eig_test - max(eps,abs(eig_test))* 10^(k-10);
            k=k+1;
          else
            find_lambda = true;
          end
        end
        neg_num_low = negative_num;
        lambda_test_low = lambda_test;    
    end

    lambda_test = eig_test;
    find_lambda = false; k = 0;
    if index == size(A,1)
        neg_zero_num_upper = size(A,1);
        lambda_test_upper = eig_test;
    else
    
        while ~find_lambda 
          if(k>10)
                error('Too many loops for finding rough upper bound')
          end
          E = mid(A) - lambda_test*mid(B);
          [L,D,P] = ldl(E);
          [negative_num, positive_num, zero_num,F] = GetInertia(D);      
          if negative_num+zero_num < index || F~=0 
            lambda_test = eig_test + max(eps,abs(eig_test)) * 10^(k-10);    
            k=k+1;
          else
            find_lambda = true;
          end
        end
        neg_zero_num_upper = negative_num+zero_num;
        lambda_test_upper = lambda_test;
    end

   %  Seek high-precision bounds. %

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Find lower bound for index-th eigenvalue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda_test = lambda_test_low;
    E = mid(A) - lambda_test*mid(B);
    [L,D,P] = ldl(mid(E));
    E = intval(A) - lambda_test*intval(B);
       
    %%%%%%%%%%%%%%  Error estimation for lower bound %%%%%%%%%%%%%%

    diff_m = P*intval(L)*intval(D)*intval(L')*P' - E;
    diff_m_mid = P*(L)*(D)*(L')* P' - mid(E);

    diff_m_inf=norm(sup(abs(diff_m_mid)),'inf');
    
    if diff_m_inf ==0
        lambda_lower = lambda_test;        
    else
        err_est = max(eps, diff_m_inf /min_eig_B);
        is_positive=0;
        k=0;
        while( k<53 && is_positive<1)
            err_est = err_est*2;

            tmpMatrix = intval(err_est)*intval(B) - diff_m;	
            is_positive = isspd( hull(tmpMatrix , tmpMatrix') );

            tmpMatrix = intval(err_est)*intval(B) + diff_m;	
            is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );

            clear tmpMatrix;
            k=k+1 ;
        end
        if( is_positive<1 )
            error('Too many loops for finding lower bound')
        end
        lambda_lower = inf( intval(lambda_test) - intval(err_est));
    end

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Find upper bound for index-th eigenvalue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    lambda_test = lambda_test_upper;
    E = mid(A) - lambda_test*mid(B);
    [L,D,P] = ldl(mid(E));
    
    E = intval(A) - lambda_test*intval(B);


    %%%%%%%%%%%%%% Error estimate for upper bound %%%%%%%%%%%%%%

    diff_m = P*intval(L)*intval(D)*intval(L')* P' - E;
    diff_m_mid = P*(L)*(D)*(L')* P' - mid(E);

    diff_m_inf = norm(sup(abs(diff_m_mid)),'inf');
    if( diff_m_inf ==0)
        lambda_upper = lambda_test;
    else
        err_est = max(eps, diff_m_inf /min_eig_B);

        is_positive=0;
        k=0;
        while( k<53 && is_positive<1)
            err_est = err_est*2;

            tmpMatrix = intval(err_est)*intval(B) - diff_m;	
            is_positive = isspd( hull(tmpMatrix , tmpMatrix'));

            tmpMatrix = intval(err_est)*intval(B) + diff_m;	
            is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );

            clear tmpMatrix;
        end
        if( is_positive<1 )
            error('Too many loops for finding upper bound')
        end

        lambda_upper = sup( intval(lambda_test) + intval(err_est) );
    end
   %----------------------------------------------------------------------
   
    eig_bound = infsup(lambda_lower, lambda_upper);

    index_range=(neg_num_low+1):1:neg_zero_num_upper;
    
    if( neg_num_low+1 < min_ind ) 
        min_ind = neg_num_low+1;
    end
    
    eig_bounds( index_range-min_ind+1,1 ) = eig_bound;
    
    index = neg_zero_num_upper+1;

end

ind_range = min_ind: (index-1); 

end

