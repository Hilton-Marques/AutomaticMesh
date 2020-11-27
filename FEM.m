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
            u = K\F;
            u = u(1:this.N,1);
        end
    end
end