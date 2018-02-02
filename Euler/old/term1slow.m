function term1 = term1slow(u,k,N,M)
% calculates term1 using an explicity convolution sum for comparison

range = [1:N,2*M-N+2:2*M];

term1 = zeros(size(u));
u = u_fullify(u,M);

k_ref = u_squishify(k,N);

for a = 1:4
    for i = 1:N
        for j = 1:N
            for l = 1:N
                current_k = squeeze(k_ref(i,j,l,:,a));
                if isequal(current_k,[0;0;0]) && ~isequal([i;j;l],[1;1;1])
                    term1(i,j,l,:,a) = [0;0;0];
                else
                    current_k
                    for i1 = range
                        for j1 = range
                            for l1 = range
                                k1 = squeeze(k(i1,j1,l1,:));
                                for i2 = range
                                    for j2 = range
                                        for l2 = range
                                            k2 = squeeze(k(i2,j2,l2,:));
                                            if k1 + k2 == current_k
                                                term1(i,j,l,:,a) = squeeze(term1(i,j,l,:,a)) - 1i*(current_k.'*squeeze(u(i1,j1,l1,:)))*squeeze(u(i2,j2,l2,:));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

