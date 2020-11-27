classdef Singular < handle
    properties
        pt
        tri
        idxPoincare
        pe1
        pe2
        pe3
        ve1
        ve2
        ve3
        t1
        t2
        t3
        u
        triSing
        color
        streams = cell.empty;
        ptBoundEdge
        ptBound
    end
    methods
        function this = Singular()
            if nargin > 0
            end
        end
        function initialize(this)
            edge1 = [this.tri.pts(1,:);this.tri.pts(2,:)];
            edge2 = [this.tri.pts(2,:);this.tri.pts(3,:)];
            edge3 = [this.tri.pts(3,:);this.tri.pts(1,:)];
            v1 = [this.u(this.tri.nodes(1),:);this.u(this.tri.nodes(2),:)];
            v2 = [this.u(this.tri.nodes(2),:);this.u(this.tri.nodes(3),:)];
            v3 = [this.u(this.tri.nodes(3),:);this.u(this.tri.nodes(1),:)];
            [this.pe1,this.ve1,this.t1] = this.intersec(edge1,v1);
            [this.pe2,this.ve2,this.t2] = this.intersec(edge2,v2);
            [this.pe3,this.ve3,this.t3] = this.intersec(edge3,v3);
        end
        function integrate(this)
            ptBound = [];
            ptBoundEdge = [];
            beginPts = {this.pe1,this.pe2,this.pe3};
            beginVs = {this.ve1,this.ve2,this.ve3};
            beginT = {this.t1,this.t2,this.t3};
            for i = 1:3
                %exportgraphics(gca,'myplot10.png','Resolution',1000); 
                pei = beginPts{i};
                vei = beginVs{i};
                tei = beginT{i};
                if ~isempty(pei)
                    for j = 1:size(pei,1)
                        stream = [];
                        pi = pei(j,:);
                        di = vei(j,:);
                        ti = tei(j);
                        triBef = this.tri;
                        edge = triBef.nodes([i,mod(i,3)+1]);
                        triNext = triBef.getAdj(i);
                        xbef = this.pt;
                        xnext = pi;
                        while (~ismember(triNext.nodes,this.triSing,'rows'))
                            %triNext.show(this.color);
                            stream(end+1,1:2) = xbef;
                            stream(end+1,1:2) = xnext;
                            line([xbef(1),xnext(1)],[xbef(2),xnext(2)],'Color','r');
                            cross1 = this.createCross(this.u(triNext.nodes(1),:));
                            cross2 = this.createCross(this.u(triNext.nodes(2),:));
                            cross3 = this.createCross(this.u(triNext.nodes(3),:));
                            cross = {cross1,cross2,cross3};
                            [~,loc] = ismember(edge,triNext.nodes);
                            if (loc(2) - loc(1) == 1)
                                edgeNexti = loc;
                                edgeNextj = [mod(edgeNexti(1),3)+1,mod(edgeNexti(2),3)+1];
                                edgeNextk = [mod(edgeNextj(1),3)+1,mod(edgeNextj(2),3)+1];
                                tnext = ti;
                            else
                                edgeNexti = [loc(2),loc(1)];
                                edgeNextj = [mod(edgeNexti(1),3)+1,mod(edgeNexti(2),3)+1];
                                edgeNextk = [mod(edgeNextj(1),3)+1,mod(edgeNextj(2),3)+1];
                                tnext = 1 - ti;
                            end
                            crossi = this.findCross(tnext,triNext.nodes,edgeNexti);
                            fxi = this.findNearestVector(crossi,di);
                            [ti,edge] = this.intercept(pi,fxi,edgeNextj,edgeNextk,triNext);
                            crossij = cross(edge);
                            fyi = this.findNearestVector(crossij{1},fxi);
                            fyj = this.findNearestVector(crossij{2},fxi);
                            fxi_ = (1-ti)*fyi + ti*fyj;
                            di = (fxi + fxi_)/2;
                            [ti,edge] = this.intercept(pi,di,edgeNextj,edgeNextk,triNext);
                            xbef = pi;
                            pi = (1 - ti)*triNext.pts(edge(1),:) + ti*triNext.pts(edge(2),:);
                            triBef = triNext;
                            triNext = triNext.getAdj(edge(1));
                            if isempty(triNext)
                                xnext = pi;
                                stream(end+1,1:2) = xbef;
                                stream(end+1,1:2) = xnext;
                                line([xbef(1),xnext(1)],[xbef(2),xnext(2)],'Color','black');
                                ptBoundEdge(end+1,1:2) = triBef.nodes(edge);
                                ptBound(end+1,1:2) = xnext;
                                %exportgraphics(gca,'myplot10.png','Resolution',1000); 
                                break;
                            end
                            edge = triBef.nodes(edge);
                            xnext = pi;
                        end
                        this.streams(end+1) = {stream};
                    end
                end
            end
            this.ptBoundEdge = ptBoundEdge;
            this.ptBound = ptBound;
        end
        function [pei,vei,tei] = intersec(this,edge,v)
            n = 20;
            pei = [];
            vei = [];
            tei = [];
            pi = edge(1,:);
            pj = edge(2,:);
            vi = v(1,:);
            vj = v(2,:);
            t = linspace(0,1,n+1);
            pk = (1 - t(1))*pi + t(1)*pj;
            u = (pk - this.pt)/norm(pk - this.pt);
            vk = (1 - t(1))*vi + t(1)*vj;
            crossk = this.createCross(vk);
            cross1bef = crossk(1,:);
            cross2bef = crossk(2,:);
            r1bef = det([cross1bef;u]);
            r2bef = det([cross2bef;u]);
            for i = 2:n+1
                pk = (1 - t(i))*pi + t(i)*pj;
                u = (pk - this.pt)/norm((pk - this.pt));
                vk = (1 - t(i))*vi + t(i)*vj;
                [r1next,r2next,cross1next,cross2next] = ...
                    this.verifySign(vk,u,cross1bef,cross2bef);
                bool1 = sign(r1bef) ~= sign(r1next);
                bool2 = sign(r2bef) ~= sign(r2next);
                if (bool1 || bool2)
                    if bool1
                        rnext = r1next;
                        crossnext = cross1next;
                    else
                        rnext = r2next;
                        crossnext = cross2next;
                    end
                    tbef = t(i-1);
                    tnext = t(i);
                    % Bissection method
                    for j = 1:10
                        rbef = rnext;
                        crossbef = crossnext;
                        newt = (tbef + tnext)/2;
                        pk = (1 - newt)*pi + newt*pj;
                        u = pk - this.pt;
                        vk = (1 - newt)*vi + newt*vj;
                        cross = this.createCross(vk);
                        crossnext = this.findNearestVector(cross,crossbef);
                        rnext = det([crossnext;u]);
                        if sign(rnext) ~= sign(rbef)
                            tbef = newt;
                        else
                            tnext = newt;
                        end
                    end
                    tk = (tnext + tbef)/2;
                    pk = (1 - tk)*pi + tk*pj;
                    u = (pk - this.pt);
                    vk = (1 - tk)*vi + tk*vj;
                    cross = this.createCross(vk);
                    vk = this.findNearestVector(cross,crossbef);
                    if (dot(u,vk) > 0 )
                        d = vk;
                    else
                        d = -vk;
                    end
                    pei(end+1,1:2) = pk;
                    tei(end+1) = tk;
                    vei(end+1,1:2) = d;
                end
                r1bef = r1next;
                r2bef = r2next;
            end
        end
        function out = findCross(this,tNext,triNextNodes,edgeNext)
            vi = triNextNodes(edgeNext(1));
            vj = triNextNodes(edgeNext(2));
            vRpr = (1 - tNext)*this.u(vi,:) + tNext*this.u(vj,:);
            out = this.createCross(vRpr);
        end
        function out = findNearestVector(~,cross,d)
            dots = zeros(1,4);
            for i = 1:4
                vk = cross(i,:);
                dots(i) = dot(vk,d);
            end
            [~,i] = max(dots);
            out = cross(i,:);
        end
        function crossVF = createCross(~,v)
            teta = mod(mod(atan2(v(:,2),v(:,1)),2*pi),2*pi);
            teta = teta/4;
            v0 = [cos(teta),sin(teta)];
            v1 = [cos(teta + pi/2),sin(teta + pi/2)];
            v2 = [cos(teta + 2*pi/2),sin(teta + 2*pi/2)];
            v3 = [cos(teta + 3*pi/2),sin(teta + 3*pi/2)];
            crossVF = [v0;v1;v2;v3];
        end
        function [r1next,r2next,cross1,cross2] = verifySign(this,vk,u,cross1bef,cross2bef)
            crossk = this.createCross(vk);
            cross1 = this.findNearestVector(crossk,cross1bef);
            cross2 = this.findNearestVector(crossk,cross2bef);
            r1next = det([cross1;u]);
            r2next = det([cross2;u]);
        end
        function [t,edge] = intercept(this,pi,vk,edge1,edge2,triNext)
            t = [];
            edge = [];
            L1 = norm(triNext.pts(edge1(2),:) - triNext.pts(edge1(1),:));
            L2 = norm(triNext.pts(edge2(2),:) - triNext.pts(edge2(1),:));
            h = L1 + L2;
            edges = {edge1,edge2};
            C = pi;
            D = pi + h*vk;
            for i = 1:2
                edge = edges{i};
                A = triNext.pts(edge(1),:) ;
                B = triNext.pts(edge(2),:);
                if this.orient(C,D,A) > 0 && this.orient(C,D,B) > 0
                    continue
                elseif this.orient(C,D,A) < 0 && this.orient(C,D,B) < 0
                    continue
                elseif this.orient(A,B,C) > 0 && this.orient(A,B,D) > 0
                    continue
                elseif this.orient(A,B,C) < 0 && this.orient(A,B,D) < 0
                end
                t = this.orient(C,D,A)/(this.orient(C,D,A) - this.orient(C,D,B));
                return
            end
        end
        function out = orient(~,A,B,C)
            out = det([B-A;C-A]);
        end
    end
    
end