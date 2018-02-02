function term2 = term2slow(u,k,N,M)
% calculates term2 using an explicit convolution sum for comparison

range = [1:N,2*M-N+2:2*M];

term2 = zeros(size(u));
u = u_fullify(u,M);

k_ref = u_squishify(k,N);

for a = 1:4
    for i = 1:N
        for j = 1:N
            for l = 1:N
                current_k = squeeze(k_ref(i,j,l,:,a));
                if isequal(current_k,[0;0;0])
                    term2(i,j,l,:,a) = [0;0;0];
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
                                                term2(i,j,l,:,a) = squeeze(term2(i,j,l,:,a)) + 1i*current_k/norm(current_k)^2*(current_k.'*squeeze(u(i1,j1,l1,:)))*(current_k.'*squeeze(u(i2,j2,l2,:)));
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

