classdef FEM < handle
    properties
        meshTri
        pts
        N
        T
        iBd
        nB
    end
    methods
        function this = FEM(meshTri,pts,nB)
            this.meshTri = meshTri;
            this.pts  = pts;
            this.iBd = 1:nB;
            this.N = size(pts,1);
            this.T = size(meshTri,1);
            this.nB = nB;
        end
        function u = Solver(this,bd)
            K=zeros(this.N,this.N);
            F=zeros(this.N,1);
            Pe = zeros(3,3);
            for i = 1: this.T
                nodes = this.meshTri(i,:);
                Pe = [ones(3,1),this.pts(nodes,:)];
                Area = abs(det(Pe))/2;
                C = inv(Pe);
                grad = C(2:3,:); Ke = Area*grad'*grad;
                K(nodes,nodes)=K(nodes,nodes)+Ke;
            end
            A = eye(this.nB);
            K(end+1:end+this.nB,1:this.nB) = A;
            K(1:this.nB,end+1:end+this.nB) = A';
            F = [F;bd];
            %u = this.cg(K,F);
            u = K\F;
            u = u(1:this.N,1);
        end
        function x = cg(~,A,b)
            max_it = 10000;
            tol = 10^-5;
            x = b;
            %  -- Iterative template routine --
            %     Univ. of Tennessee and Oak Ridge National Laboratory
            %     October 1, 1993
            %     Details of this algorithm are described in "Templates for the
            %     Solution of Linear Systems: Building Blocks for Iterative
            %     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
            %     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
            %     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
            %
            %  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
            %
            % cg.m solves the symmetric positive definite linear system Ax=b
            % using the Conjugate Gradient method with preconditioning.
            %
            % input   A        REAL symmetric positive definite matrix
            %         x        REAL initial guess vector
            %         b        REAL right hand side vector
            %         M        REAL preconditioner matrix
            %         max_it   INTEGER maximum number of iterations
            %         tol      REAL error tolerance
            %
            % output  x        REAL solution vector
            %         error    REAL error norm
            %         iter     INTEGER number of iterations performed
            %         flag     INTEGER: 0 = solution found to tolerance
            %                           1 = no convergence given max_it
            
            flag = 0;                                 % initialization
            iter = 0;
            
            bnrm2 = norm( b );
            if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
            
            r = b - A*x;
            error = norm( r ) / bnrm2;
    
            
            for iter = 1:max_it                       % begin iteration
                
                z  = r;
                rho = (r'*z);
                
                if ( iter > 1 )                       % direction vector
                    beta = rho / rho_1;
                    p = z + beta*p;
                else
                    p = z;
                end
                
                q = A*p;
                alpha = rho / (p'*q );
                x = x + alpha * p;                    % update approximation vector
                
                r = r - alpha*q;                      % compute residual
                error = norm( r ) / bnrm2;            % check convergence
                if ( error <= tol )
                    break
                end
                
                rho_1 = rho;
                
            end
            
            if ( error > tol ) 
                flag = 1; 
            end         % no convergence
        end
            % END cg.m
        end
    end