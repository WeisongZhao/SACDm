function f=imreadstack(imname)
info = imfinfo(imname);
num_images = numel(info);
f=zeros(info(1).Height,info(1).Width,num_images);

for k = 1:num_images
    f(:,:,k) =single(imread(imname, k));
end