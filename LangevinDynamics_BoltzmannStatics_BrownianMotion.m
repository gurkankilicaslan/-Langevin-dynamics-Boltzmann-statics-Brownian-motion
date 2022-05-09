tic
clear; clc; clf;
T = 1000;
L = [50 50]; 
N = 1000; 
D = 1/4; 
r(1:N,1) = rand(N,1)*(L(1,1)-D/2); 
r(1:N,2) = rand(N,1)*(L(1,2)-D/2);
v = randn(N,2);
dt = 0.1;

mp = 1; 
mb = 10;

vb = [0 0];
rb = [25 25];
pb(T,2) = 0;
Db = 15; 

nv = find((r(:,1) - rb(1,1)).^2 + (r(:,2) - rb(1,2)).^2 < ((Db+D)/2)^2);
for n = 1 : size(nv)
    th = rand*2*pi;
    r(nv(n),:) = r(nv(n),:) + (Db+D)*[cos(th) sin(th)]; 
end


for t = 1 : T
    summ = 0;
    pb(t,:) = rb;
    for a = 1:200
            
        

        plot(r(:,1),r(:,2),'.r','MarkerSize',20); 
        hold on;
        plot(pb(1:t,1),pb(1:t,2),"b");
        rectangle('Position',[rb(1,1)-Db/2 rb(1,2)-Db/2 Db Db],'Curvature',[1 1]); % plot particles
    
        rectangle('Position',[0 0 L(1,1) L(1,2)]);
        hold off;
        axis equal;
        axis([0 L(1,1) 0 L(1,2)]);
        drawnow;

        
        [sx,nx] = sort(r(:,1)); 
        for kn = 1 : N-1
            n = nx(kn); 
            for km = kn+1:N
                m = nx(km);
                if (sx(km) - sx(kn) > D) 
                    break
                else
                    rp = r(n,:)-r(m,:); 
                    nrp = norm(rp);
                    if (nrp < D)
                        rv = v(n,:)-v(m,:); 
                        if (rv*rp'<0)
                            v(n,:) = v(n,:) - (rv*rp')*rp/nrp^2; 
                            v(m,:) = v(m,:) + (rv*rp')*rp/nrp^2; 
                        end
                    end
                end
            end
        end

        for n = 1 : 2 
            nv = (v(:,n) > 0).*(r(:,n) > L(1,n)-D/2); v(nv==1,n) = -v(nv==1,n);
            nv = (v(:,n) < 0).*(r(:,n) < D/2); v(nv==1,n) = -v(nv==1,n);
        end

        nv = find((r(:,1) - rb(1,1)).^2 + (r(:,2) - rb(1,2)).^2 < ((Db+D)/2)^2);
        for n = 1:size(nv,1)
            m = nv(n,1);
            rp = (r(m,:)-rb); nrp = norm(rp); rv = v(m,:) - vb;
            if (rv*rp' < 0)
                v(m,:) = v(m,:) - 2*mb/(mp + mb)*(rv*rp')*rp/nrp^2;
                vb = vb + 2*mp/(mp + mb)*(rv*rp')*rp/nrp^2;
            end
        end

        
        r = r + v * dt;
        rb = rb + vb * dt; 
        
        
        
        mp*sum(v(:,1).^2 + v(:,2).^2) + mb*(vb(1,1)^2 + vb(1,2)^2);
        summ = summ + (abs(rb(1)-25))^2 + (abs(rb(2)-25))^2;
    end
    average = summ/200;
    result(t)=average;
end
plot((1:T),result)
toc