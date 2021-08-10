function Background = background_estimation(imgs,th,dlevel,wavename,iter)
if nargin < 2 || isempty(th)
    th=1;
end
if nargin < 3 || isempty(dlevel)
    dlevel=7;
end
if nargin < 4 || isempty(wavename)
    wavename='db6';
end
if nargin < 5 || isempty(iter)
    iter=3;
end
[x,y,~]=size(imgs);
if x<y
    imgs=padarray(imgs,[max(x,y)-size(imgs,1)...
        ,max(x,y)-size(imgs,2),0],'post');
end
Background = zeros(size(imgs),'single');
for frames = 1: size(imgs,3)
    initial = imgs(:,:,frames);
    res = initial;
    for ii = 1:iter
        [m,n] = wavedec2(res,dlevel,wavename);
        vec = zeros(size(m));
        vec(1:n(1)*n(1)*1) = m(1:n(1)*n(1)*1);
        Biter =  waverec2(vec,n,wavename);
        if th > 0
            eps = sqrt(abs(res))/2;
            ind = initial>(Biter+eps);
            res(ind) = Biter(ind)+eps(ind);
            [m,n] = wavedec2(res,dlevel,wavename);
            vec = zeros(size(m));
            vec(1:n(1)*n(1)*1) = m(1:n(1)*n(1)*1);
            Biter =  waverec2(vec,n,wavename);
        end
    end
    Background(:,:,frames) = Biter;
end
Background=Background(1:x,1:y,:);