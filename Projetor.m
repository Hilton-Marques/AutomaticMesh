classdef Projetor < handle
    properties
        Ci = cell.empty;
        meshx ;
        meshy ;
        nu;
        nv;
    end
    methods
        function this = Projetor(Ci)
            C1 = Ci{1};
            C2 = Ci{2};
            C3 = Ci{3};
            C4 = Ci{4};
            nC1 = size(C1,1);
            nC2 = size(C2,1);
            nC3 = size(C3,1);
            nC4 = size(C4,1);
            Nmax = min([nC1,nC2,nC3,nC4]);
            Nmax = 10;
            Ci{1} = this.correction(C1,Nmax);
            Ci{2} = this.correction(C2,Nmax);
            Ci{3} = this.correction(C3,Nmax);
            Ci{4} = this.correction(C4,Nmax);
            %             if nC1 > nC2
            %                 Ci{2} = this.correction(C2,C1);
            %             elseif nC2 > nC1
            %                 Ci{1} = this.correction(C1,C2);
            %             end
            %             if nC3 > nC4
            %                 Ci{4} = this.correction(C4,C3);
            %             elseif nC4 > nC3
            %                 Ci{3} = this.correction(C3,C4);
            %             end
            this.Ci = Ci;
            this.start()
        end
        function start(this)
            Ci = this.Ci;
            nu = size(Ci{1},1);
            nv = size(Ci{3},1);
            u = linspace(0,1,nu);
            v = linspace(0,1,nv);
%             % slow form
%             mmx = zeros(nv,nu);
%             for j = 1: nu
%                 for i = 1:nv
%                     mmx(i,j) = (1 - v(i))*Ci{1}(j,1) + ...
%                         v(i)*Ci{2}(j,1) + (1-u(j))*Ci{3}(i,1) + ...
%                         u(j)*Ci{4}(i,1) - (1 - u(j))*( 1- v(i))*Ci{1}(1,1) - ...
%                         u(j)*(1 - v(i))*Ci{1}(end,1) - u(j)*v(i)*Ci{2}(end,1) ...
%                          - (1 - u(j))*v(i)*Ci{2}(1,1);
%                 end
%             end
            pvx = (1-v)'.*Ci{1}(:,1)' + v'.*Ci{2}(:,1)';
            pvy = (1-v)'.*Ci{1}(:,2)' + v'.*Ci{2}(:,2)';
            pux = (1-u).*Ci{3}(:,1) + u.*Ci{4}(:,1);
            puy = (1-u).*Ci{3}(:,2) + u.*Ci{4}(:,2);
            ppx = (1 - u) .* (1 - v)' * Ci{1}(1,1) + ...
                u.*(1 - v)'*Ci{1}(end,1) + u.*v'*Ci{2}(end,1) + ...
                (1-u).* v' * Ci{2}(1,1);
            ppy = (1 - u) .* (1 - v)' * Ci{1}(1,2) + ...
                u.*(1 - v)'*Ci{1}(end,2) + u.*v'*Ci{2}(end,2) + ...
                (1-u).* v' * Ci{2}(1,2);
            
            this.meshx = pvx + pux  - ppx ;
            this.meshy = pvy + puy - ppy ;
            this.nu = nu;
            this.nv = nv;
        end
        function plot(this)
            meshx = this.meshx;
            meshy = this.meshy;
            nv = this.nv;
            nu = this.nu;
            for i = 1:nv-1
                for j = 1: nu-1
                    p1 = [meshx(i,j),meshy(i,j)];
                    p2 = [meshx(i,j+1),meshy(i,j+1)];
                    p3 = [meshx(i+1,j+1),meshy(i+1,j+1)];
                    p4 = [meshx(i+1,j),meshy(i+1,j)];
                    px = [p1(1),p2(1),p3(1),p4(1)];
                    py = [p1(2),p2(2),p3(2),p4(2)];
                    this.drawRect(p1,p2,p3,p4);
                end
            end
        end
        function newCMenor = correction(this,curve,nc_maior)
            L = this.totalLength(curve);
            newCMenor = zeros(nc_maior,2);
            newCMenor(1,:) = curve(1,:);
            newCMenor(end,:) = curve(end,:);
            Lu = L/(nc_maior-1);
            for i = 1: nc_maior-2
                target = i*Lu;
                [id,h] = this.accumulate(curve,target);
                a = curve(id - 1,:);
                b = curve(id,:);
                p = a + h*(b-a);
                newCMenor(i+1,:) = p;
            end
        end
        function drawRect(~,p1,p2,p3,p4)
            line([p1(1),p2(1)],[p1(2),p2(2)]);
            line([p2(1),p3(1)],[p2(2),p3(2)]);
            line([p3(1),p4(1)],[p3(2),p4(2)]);
            line([p4(1),p1(1)],[p4(2),p1(2)]);
        end
        function out = totalLength(~,curve)
            N = size(curve,1);
            out = 0;
            p1 = curve(1,:);
            for i = 2:N
                p2 = curve(i,:);
                Li = norm(p2 - p1);
                out = Li + out;
                p1 = p2;
            end
        end
        function [id,h] = accumulate(~,curve,target)
            N = size(curve,1);
            L = 0;
            p1 = curve(1,:);
            for i = 2:N
                p2 = curve(i,:);
                Li = norm(p2 - p1);
                L = Li + L;
                p1 = p2;
                if L >= target
                    id = i;
                    h = (target - (L - Li))/Li;
                    return
                end
            end
        end
    end
end