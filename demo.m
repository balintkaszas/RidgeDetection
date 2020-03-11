domain = [0,6;-34,-28];
resolution = [400,400];

dataset = load('data/ftle.mat');
ftles = dataset.ftle_;
C = imgaussfilt(ftles,3); %%Gaussian filter with stdev = 3;

[x,y] = detectRidge(C, resolution, domain);
plot(x,y, '.')
xlim(domain(1,:))
ylim(domain(2,:))