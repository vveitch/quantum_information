function [ plot ] = plot_hack( rep_mat)
%PLOT_HACK plots the wig rep of a density matrix
%   a hack to get around the irritating interpolation thing of pcolor

inflat = 10;

hack_rep=zeros(inflat*size(rep_mat,1));
for i=0:size(rep_mat,1)-1
    for j=0:size(rep_mat,2)-1
    hack_rep(inflat*i+1:inflat*(i+1), inflat*j+1:inflat*(j+1))=rep_mat(i+1,j+1);
    end
end
hack_rep
plot=pcolor(hack_rep);

end

