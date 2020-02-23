domain = [0,2;0,1];
resolution = [100,100];
difference = 3./1000;
initialPosition = initialize_ic_grid(resolution, domain, 2);
coords = reshape(initialPosition, [100,100,2]);
coords = coords(3:end-3, 3:end-3,:);
xi = coords(:,:,1);
yi = coords(:,:,2);
asd = load('doublegyre.mat');
ftles = asd.a;
ftles = ftles(3:end-3, 3:end-3);
FF = repmat(ftles,1,1,3);
FF=smooth3(FF,'gaussian',3);
smoothFT = FF(:,:,1);



% % surf(xi,yi,FF(:,:,1));shading interp; axis equal;axis tight;colorbar;
% % view([0 0 1]); axis equal; axis tight; shading interp;camlight

[fx,fy] = gradient(smoothFT, difference);
%imagesc(domain(:,1), domain(:,2), fx);
[fxx, fxy] = gradient(fx, difference);
[fyx, fyy] = gradient(fy, difference);

eigenvalsBIG = 0.5*(fxx + fyy + sqrt(fxx.^2 + 4. * fxy.^2 - 2.*fxx.*fyy + fyy.^2));
eigenvalsSmall = 0.5*(fxx + fyy - sqrt(fxx.^2 + 4. * fxy.^2 - 2.*fxx.*fyy + fyy.^2));
eigenvectorx = fyy - eigenvalsSmall;
eigenvectory = fyx;
magnitudes = sqrt(fx.^2 + fy.^2);
mask2 = eigenvalsBIG<0;
ix = find((imregionalmax(smoothFT).*mask2));

% surf(xi,yi,FF(:,:,1));shading interp; axis equal;axis tight;colorbar;
% view([0 0 1]); axis equal; axis tight; shading interp;camlight
hold on;
imagesc(domain(1,:), domain(2,:), smoothFT);
plot(xi(ix),yi(ix), 'r*','MarkerSize',20)


%imagesc(domain(1,:), domain(2,:), eigenvalsSmall);
mask = magnitudes > mean(magnitudes,  'all')*2;
%magnitudes = sqrt(fx.^2 + fy.^2);
%disp(min(magnitudes, [], 'all'));
% surf(xi,yi,magnitudes);shading interp; axis equal;axis tight;colorbar;
% view([0 0 1]); axis equal; axis tight; shading interp;camlight
%histogram(reshape(magnitudes, [1000000,1]));

% 
 %mask = magnitudes > 1e3;
%imagesc(domain(1,:), domain(2,:), mask);



sum(mask, 'all')
maskedx = xi(mask);
maskedy = yi(mask);


xi = xi';
yi = yi';
fx = fx';
fy = fy';
fxInterp = griddedInterpolant(xi,yi,fx, 'spline');
fyInterp = griddedInterpolant(xi,yi,fy, 'spline');
eigenvalueSmallInterp = griddedInterpolant(xi,yi,eigenvalsSmall, 'spline');
eigenvectorxSmallInterp = griddedInterpolant(xi,yi,eigenvectorx, 'spline');
eigenvectorySmallInterp = griddedInterpolant(xi,yi,eigenvectory, 'spline');



%Opt    = odeset('Events', @event);

% 
% %[t,sol] = ode45(@(t,e) odefun(t,e, @(x1, x2)fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2) ), [0, 0.1*pi], [maskedx(10), maskedy(10)], optionsode);
% %plot(sol(:,1), sol(:,2), '.', 'color', 'black','LineWidth', 4)
hold on;
N = sum(mask, 'all');
ridgelocations = zeros(N,2);
for i = 1:N
    disp(i/N);

    eigenvec = [eigenvectorxSmallInterp(maskedx(i),maskedy(i)), eigenvectorySmallInterp(maskedx(i),maskedy(i))];
    gradie = [fxInterp(maskedx(i),maskedy(i)), fyInterp(maskedx(i),maskedy(i))];
    prevangle = dot(gradie,eigenvec)/(norm(gradie)*norm(eigenvec));
    eventfonction = @(t, y) eventss(t, y, @(x1, x2) eigenvalueSmallInterp(x1, x2), @(x1, x2) eigenvectorxSmallInterp(x1, x2), @(x1, x2) eigenvectorySmallInterp(x1, x2), @(x1, x2) fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2), prevangle);
    optionsode=odeset('Events',eventfonction);
    if((eigenvalueSmallInterp(maskedx(i),maskedy(i)) > 0. ))
    [~,sol] = ode45(@(t,e) odefun(t,e, @(x1, x2)fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2) ), [0, 10], [maskedx(i), maskedy(i)], optionsode);
     plot(sol(:,1), sol(:,2), '.', 'color', 'black','LineWidth', 4)
        if(~isnan(sol(end)))
            ridgelocations(i,:) = sol(end, :);
        end
        
    end
    
end
%eigenvalueSmallInterp(3, 4);
function [value, isterminal, direction] = eventss(T, X, eigenvalueSmallInterp, eigenvectorxSmallInterp, eigenvectorySmallInterp, fxInterp, fyInterp, prevangle)
x = X(1);
y = X(2);

eigenvec = [eigenvectorxSmallInterp(x,y), eigenvectorySmallInterp(x,y)];
gradie = [fxInterp(x,y), fyInterp(x,y)];
angle = dot(gradie,eigenvec)/(norm(gradie)*norm(eigenvec));

value      = (eigenvalueSmallInterp(x,y) > 0. ) && ( abs(angle-prevangle) > 1e-5) && ~isnan(x) && isInRegion(X);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function dy = odefun(t, X, fx, fy)
x = X(1);
y = X(2);

dy(1) = fx(x, y);

dy(2) = fy(x, y);
dy = dy';
end

function isinReg = isInRegion(x)
isinReg = 1;
if(x(1)<0 || x(1)>2 || x(2) > 1 || x(1)<0)
    isinReg = 0;
end
end