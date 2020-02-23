domain = [-1.5, 1.5; -1.5, 1.5];
resolution = [200,200];
difference = 3./200;
initialPosition = initialize_ic_grid(resolution, domain, 2);

%domain = [-1.5,-0.9;-1, -0.4];
domain = [-1.5, 1.5;-1.5, 1.5];
coords = reshape(initialPosition, [200,200,2]);
xi = coords(:,:,1);
yi = coords(:,:,2);
asd = load('Duff_cons_200by200_0_4pi.mat');
ftles = asd.ftle;
FF = repmat(ftles,1,1,3);
FF=smooth3(FF,'gaussian',3);
smoothFT = FF(:,:,1);



[fx,fy] = gradient(smoothFT, difference);
%imagesc(domain(:,1), domain(:,2), fx);
[fxx, fxy] = gradient(fx, difference);
[fyx, fyy] = gradient(fy, difference);

eigenvalsBIG = 0.5*(fxx + fyy + sqrt(fxx.^2 + 4. * fxy.^2 - 2.*fxx.*fyy + fyy.^2));
eigenvalsSmall = 0.5*(fxx + fyy - sqrt(fxx.^2 + 4. * fxy.^2 - 2.*fxx.*fyy + fyy.^2));
eigenvectorx = fyy - eigenvalsSmall;
eigenvectory = fyx;
magnitudes = sqrt(fx.^2 + fy.^2);


surf(xi,yi,magnitudes);shading interp; axis equal;axis tight;colorbar;
view([0 0 1]); axis equal; axis tight; shading interp;camlight


%imagesc(domain(1,:), domain(2,:), eigenvalsSmall);
mask = magnitudes > max(magnitudes, [], 'all')*0.9;
%magnitudes = sqrt(fx.^2 + fy.^2);
%disp(min(magnitudes, [], 'all'));
%surf(xi,yi,magnitudes);shading interp; axis equal;axis tight;colorbar;
%view([0 0 1]); axis equal; axis tight; shading interp;camlight
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


hold on;
%imagesc(domain(1,:), domain(2,:), ftles);
%Opt    = odeset('Events', @event);
% eventfonction = @(t, y) eventss(t, y, @(x1, x2) eigenvalueSmallInterp(x1, x2), @(x1, x2) eigenvectorxSmallInterp(x1, x2), @(x1, x2) eigenvectorySmallInterp(x1, x2), @(x1, x2) fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2));
% optionsode=odeset('Events',eventfonction);
% 
% %[t,sol] = ode45(@(t,e) odefun(t,e, @(x1, x2)fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2) ), [0, 0.1*pi], [maskedx(10), maskedy(10)], optionsode);
%plot(sol(:,1), sol(:,2), '.', 'color', 'black','LineWidth', 4)
% 
% for i = 5:17
%     [~,sol] = ode45(@(t,e) odefun(t,e, @(x1, x2)fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2) ), [0, pi], [maskedx(i), maskedy(i)], optionsode);
%     %plot(maskedx, maskedy, 'o');
%     plot(sol(:,1), sol(:,2), '.', 'color', 'black','LineWidth', 4)
%     %if(norm(sol(end)) > 1e10)
%         %disp(norm(sol(end)))
%     %end
% end
%eigenvalueSmallInterp(3, 4);
function [value, isterminal, direction] = eventss(T, X, eigenvalueSmallInterp, eigenvectorxSmallInterp, eigenvectorySmallInterp, fxInterp, fyInterp)
x = X(1);
y = X(2);

eigenvec = [eigenvectorxSmallInterp(x,y), eigenvectorySmallInterp(x,y)];
gradie = [fxInterp(x,y), fyInterp(x,y)];
disp(eigenvalueSmallInterp(x,y));
disp(dot(gradie,eigenvec)/(norm(gradie)*norm(eigenvec)));

value      = (eigenvalueSmallInterp(x,y) > 0. ) && ( abs(dot(gradie,eigenvec)/(norm(gradie)*norm(eigenvec))-1.) > 1e-2);
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