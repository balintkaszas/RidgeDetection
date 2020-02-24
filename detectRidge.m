function ridgeMask = detectRidge(scalarField, resolution, domain)
%% Ridge detection, by computing the gradient flow of the scalar field.


%computing the gradient:

%grid spacing in either direction
differenceX = diff(domain(1,:))/resolution(1);
differenceY = diff(domain(2,:))/resolution(2);


%generating the uniform grid
initialPosition = initialize_ic_grid(resolution, domain, 2);
coords = reshape(initialPosition, [resolution(1), resolution(2), 2]);
%coords = coords(3:end-3, 3:end-3,:);
xi = coords(:,:,1);
yi = coords(:,:,2);
%asd = load('doublegyre.mat');
%ftles = asd.a;
%ftles = ftles(3:end-3, 3:end-3);
FF = repmat(scalarField,1,1,3); %%apply initial smoothing to the scalarField
FF=smooth3(FF,'gaussian',3);
smoothScalarField = FF(:,:,1);
 

%calculate the gradient and smooth each component
[fx,fy] = gradient(smoothScalarField, differenceX, differenceY);
fx = repmat(fx,1,1,3);
fx=smooth3(fx,'gaussian',3);
smoothfx = fx(:,:,1);
fy = repmat(fy,1,1,3);
fy=smooth3(fy,'gaussian',3);
smoothfy = fy(:,:,1);


% calculate Hessian
[fxx, fxy] = gradient(smoothfx, differenceX, differenceY);
[fyx, fyy] = gradient(smoothfy, differenceX, differenceY);

% smooth components of the Hessian
fxx = repmat(fxx,1,1,3);
fxx=smooth3(fxx,'gaussian',3);
smoothfxx = fxx(:,:,1);

fyy = repmat(fyy,1,1,3);
fyy=smooth3(fyy,'gaussian',3);
smoothfyy = fyy(:,:,1);

fxy = repmat(fxy,1,1,3);
fxy=smooth3(fxy,'gaussian',3);
smoothfxy = fxy(:,:,1);

fyx = repmat(fyx,1,1,3);
fyx=smooth3(fyx,'gaussian',3);
smoothfyx = fyx(:,:,1);


% calculate eigenvalues and eigenvectors of the hessian explicitly
eigenvalsBIG = 0.5*(smoothfxx + smoothfyy + sqrt(smoothfxx.^2 + 4. * smoothfxy.^2 - 2.*smoothfxx.*smoothfyy + smoothfyy.^2));
eigenvalsSmall = 0.5*(smoothfxx + smoothfyy - sqrt(smoothfxx.^2 + 4. * smoothfxy.^2 - 2.*smoothfxx.*smoothfyy + smoothfyy.^2));


eigenvectorx = smoothfyy - eigenvalsSmall;
eigenvectory = smoothfyx;

% Step1: find points which have a high gradient
magnitudes = sqrt(smoothfx.^2 + smoothfy.^2);
mask = magnitudes > mean(magnitudes,  'all');
mask(1:3,:) = 0;
mask(end-3:end,:) = 0;
mask(:,1:3) = 0;
mask(:,end-3:end) = 0;

%ix = find((imregionalmax(smoothFT)));


maskedx = xi(mask);
maskedy = yi(mask);



% interpolate each function defined on the grid:
xi = xi';
yi = yi';
smoothfx = smoothfx';
smoothfy = smoothfy'; %this is switching to ndview
%quiver(xi,yi,smoothfx,smoothfy,'r'); figure(gcf)
fxInterp = griddedInterpolant(xi,yi,smoothfx, 'spline');
fyInterp = griddedInterpolant(xi,yi,smoothfy, 'spline');
eigenvalueSmallInterp = griddedInterpolant(xi,yi,eigenvalsSmall, 'spline');

eigenvectorxSmallInterp = griddedInterpolant(xi,yi,eigenvectorx, 'spline');
eigenvectorySmallInterp = griddedInterpolant(xi,yi,eigenvectory, 'spline');

N = sum(mask, 'all'); %total number of initial points
%create empty mask:
zeroMatrix = zeros(size(initialPosition));


for i = 1:N
    disp(i/N);
    %calculate the angle between the eigenvector of the hessian and the
    %gradient 
    eigenvec = [eigenvectorxSmallInterp(maskedx(i),maskedy(i)), eigenvectorySmallInterp(maskedx(i),maskedy(i))];
    gradie = [fxInterp(maskedx(i),maskedy(i)), fyInterp(maskedx(i),maskedy(i))];
    prevangle = dot(gradie,eigenvec)/(sqrt(norm(gradie)*norm(eigenvec)));
    
    eventfonction = @(t, y) eventss(t, y, @(x1, x2) eigenvalueSmallInterp(x1, x2), @(x1, x2) eigenvectorxSmallInterp(x1, x2), @(x1, x2) eigenvectorySmallInterp(x1, x2), @(x1, x2) fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2), prevangle, domain);
    optionsode=odeset('Events',eventfonction, 'RelTol', 1e-5);
    
    
    if((eigenvalueSmallInterp(maskedx(i),maskedy(i)) > 0. ))
        [~,sol] = ode45(@(t,e) odefun(t,e, @(x1, x2)fxInterp(x1, x2), @(x1, x2) fyInterp(x1, x2) ), [0, 100], [maskedx(i), maskedy(i)], optionsode);
        if(~isnan(sol(end)) && isInRegion(sol(end,:), domain)) %if the converged solution is inside the region
            %find the gridpoint closest to it
            k = dsearchn(initialPosition, sol(end,:)); %returns linear index, at which initialPosition(k,:) ~= sol(end,:).
            zeroMatrix(k,:) = 1; %set the value of the mask to 1 
        end
        
    end
end

ridgeMask = reshape(zeroMatrix(:,1), resolution);

end

%eigenvalueSmallInterp(3, 4);
function [value, isterminal, direction] = eventss(T, X, eigenvalueSmallInterp, eigenvectorxSmallInterp, eigenvectorySmallInterp, fxInterp, fyInterp, prevangle, domain)
x = X(1);
y = X(2);

eigenvec = [eigenvectorxSmallInterp(x,y), eigenvectorySmallInterp(x,y)];
gradie = [fxInterp(x,y), fyInterp(x,y)];
angle = dot(gradie,eigenvec)/(sqrt(norm(gradie)*norm(eigenvec)));

value      = (eigenvalueSmallInterp(x,y) > 0. ) && ( abs(angle-prevangle) > 1e-5) && ~isnan(x) && isInRegion(X, domain);
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

function isinReg = isInRegion(x, domain)
isinReg = 1;
    if(x(1) < domain(1,1) || x(1) > domain(1,2) || x(2) < domain(2,1) || x(2) > domain(2,2))
        isinReg = 0;
    end
end