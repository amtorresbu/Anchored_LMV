

##

## ## Functions for reconstruction --------------------------------

##

## (L1).0 ----------------------------------------------------------

##

function Rec_L1_0(m, s, ε)
    
    ## Input: m views, s point correspondences, ε noise level
    #
    ## Output: The average log relative error of reconstructing s point correspondences using (L1).0
    #

    cams = [randn(3, 4) for _ in 1:m]
    cams_vec = vec(cams[1]) 
    for t in 2:m cams_vec=vcat(cams_vec,vec(cams[t])) end 

    v = randn(2, 2); # Parameters for the random line
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0] # Random line

    k = randn(1,s)
    PP = [line*[k[j];1] for j in 1:s] #points
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s]
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s]

    ## adding noise to points
    ee=rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s]

    error_sum = 0
    
    for j in 1:s # Reconstruct all the point correspondences
        
        qq_vec = vec(qq[j][1][1:2]) 
        for t in 2:m qq_vec = vcat(qq_vec,vec(qq[j][t][1:2])) end 

        R = real_solutions(HomotopyContinuation.solve(D_PMV, start_PMV , start_parameters = PMV_q₀,  target_parameters = [cams_vec;qq_vec],show_progress=false))

        min=10 # High starting error
        original_point=PP[j]/PP[j][4,1]
        for k in 1:length(R)
            new_point=[R[k][1:3];1]
            if min>norm(original_point-new_point)
                min=norm(original_point-new_point)
            end
        end

        error_sum = error_sum+min 
    end

    return log(10,(error_sum/s)/ε)
end

function Rec_L1_0_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #         
    ## Output: The average error of reconstructing s point correspondences using (L1).0
    #
    
    errSum = []; 

    for i in 1:iter errSum=vcat(errSum,Rec_L1_0(m, s, ε)) end

    return errSum 
end

function Rec_L1_0_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences,  ε noise level
    #        
    ## Output: The average time of reconstructing s point correspondences using (L1).0
    #
    
    time = [];
    
    for i in 1:iter time = vcat(time, @elapsed Rec_L1_0(m,s,ε)) end

    return time 
end

##

## (L1).1 ----------------------------------------------------------

##

∧(M) = M[:,1] × M[:,2]; 

function Rec_L1_1(m, s, ε)

    ## Input: m views, s point correspondences, ε noise level
    #
    ## Output: The average log relative error of reconstructing s point correspondences using (L1).1
    #

    cams = [randn(3, 4) for _ in 1:m];

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0]; # Random line

    k = randn(1, s);
    PP = [line*[k[j];1] for j in 1:s]; #points
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s];
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s];
    
    ll = [(∧(cams[i]*line)) for i in 1:m];
    ll = [ll[i]/ll[i][3] for i in 1:m];

    ## adding noise to points
    ee = rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    ee = rand(2,m)
    uu = [[(ll[i][1:2] + ε*ee[:,i]/norm(ee[:,i]) );1] for i in 1:m];

    
    bp_line = nullspace(transpose([transpose(cams[1])*uu[1] transpose(cams[2])*uu[2]])) # ON basis of the line

    AA = [transpose(Matrix(qr(cams[i]*bp_line).Q)) for i in 1:m] # Projection matrices X such that XX^T=I  
    qq_A = [[AA[i]*qq[j][i] for i in 1:m] for j in 1:s] 
    qq_A = [[qq_A[j][i]/qq_A[j][i][2] for i in 1:m] for j in 1:s]

    cams_aug = [AA[i]*cams[i]*bp_line for i in 1:m]
    cams_aug_vec = vcat(vec(cams_aug[1]))
    for t in 2:m cams_aug_vec=vcat(cams_aug_vec,vec(cams_aug[t])) end 

    error_sum = 0
    
    for j in 1:s # Reconstruct all the point correspondences in the anchored PMV

        qq_A_vec = qq_A[j][1][1]
        for t in 2:m qq_A_vec=vcat(qq_A_vec, qq_A[j][t][1]) end

        R = real_solutions(HomotopyContinuation.solve(D_M11, start_M11 , start_parameters = M11_q₀,  target_parameters = [cams_aug_vec;qq_A_vec],show_progress=false))
        
        min=10
        original_point=PP[j]/PP[j][4,1]
        for k in 1:length(R)
            new_point=bp_line*[R[k][1];1]
            new_point=new_point/new_point[4,1]
            if min>norm(original_point-new_point)
                min=norm(original_point-new_point)
            end
        end

        error_sum = error_sum+min
    end

   return log(10,(error_sum/s)/ε)
end

function Rec_L1_1_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspondences using (L1).1
    #
    
    errSum = []; 

    for i in 1:iter errSum=vcat(errSum,Rec_L1_1(m, s, ε)) end

    return errSum 
end

function Rec_L1_1_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average time of reconstructing s point correspondence using (L1).1
    #
    
    time = []; 

    for i in 1:iter time=vcat(time,@elapsed Rec_L1_1(m, s, ε)) end

    return time 
end

##

## (L1).2 ----------------------------------------------------------

##


function Rec_L1_2(m, s, ε)

    ## Input: m views, s point correspondences, ε noise level
    #
    ## Output: The average log relative error of reconstructing s point correspondences using (L1).2
    #

    cams = [randn(3, 4) for _ in 1:m];
    cams_vec = vec(cams[1]) 
    for t in 2:m cams_vec=vcat(cams_vec,vec(cams[t])) end

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0]; # Random line

    k = sort([randn(1) for i in 1:s]);
    PP = [line*[k[j];1] for j in 1:s]; # Points
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s];
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s];

    ## adding noise to points
    ee = rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    # Reconstruct the two endpoints
    qq_vec_1 = vec(qq[1][1][1:2]) 
    for t in 2:m qq_vec_1 = vcat(qq_vec_1,vec(qq[1][t][1:2])) end 

    R = real_solutions(HomotopyContinuation.solve(D_PMV, start_PMV , start_parameters = PMV_q₀,  target_parameters = [cams_vec;qq_vec_1],show_progress=false))
    
    original_point_1 = PP[1]/PP[1][4,1]
    best_point_1 = [R[1][1:3];1]
    min = norm(original_point_1-best_point_1)
    for k in 2:length(R)
        new_point = [R[k][1:3];1]
        if min > norm(original_point_1-new_point)
            min = norm(original_point_1-new_point)
            best_point_1 = new_point
        end
    end

    error_sum = min
    
    qq_vec_2 = vec(qq[s][1][1:2]) 
    for t in 2:m qq_vec_2 = vcat(qq_vec_2,vec(qq[s][t][1:2])) end 

    R = real_solutions(HomotopyContinuation.solve(D_PMV, start_PMV , start_parameters = PMV_q₀,  target_parameters = [cams_vec;qq_vec_2],show_progress=false))
    
    original_point_2 = PP[s]/PP[s][4,1]
    best_point_2 = [R[1][1:3];1]
    min = norm(original_point_2-best_point_2)
    for k in 2:length(R)
        new_point = [R[k][1:3];1]
        if min > norm(original_point_2-new_point)
            min = norm(original_point_2-new_point)
            best_point_2 = new_point
        end
    end

    error_sum = error_sum + min

    bp_line = Matrix(qr(hcat(best_point_1, best_point_2)).Q) # ON basis of the line
    
    AA = [transpose(Matrix(qr(cams[i]*bp_line).Q)) for i in 1:m] # Projection matrices X such that XX^T=I  
    qq_A = [[AA[i]*qq[j][i] for i in 1:m] for j in 1:s] 
    qq_A = [[qq_A[j][i]/qq_A[j][i][2] for i in 1:m] for j in 1:s]

    cams_aug = [AA[i]*cams[i]*bp_line for i in 1:m]
    cams_aug_vec = vcat(vec(cams_aug[1]))
    for t in 2:m cams_aug_vec=vcat(cams_aug_vec,vec(cams_aug[t])) end 
    
    for j in 2:s-1 # Reconstruct all the point correspondences in the anchored PMV

        qq_A_vec = qq_A[j][1][1]
        for t in 2:m qq_A_vec=vcat(qq_A_vec, qq_A[j][t][1]) end 

        R = real_solutions(HomotopyContinuation.solve(D_M11, start_M11 , start_parameters = M11_q₀,  target_parameters = [cams_aug_vec;qq_A_vec],show_progress=false))
        
        min=10
        original_point=PP[j]/PP[j][4,1]
        for k in 1:length(R)
            new_point=bp_line*[R[k][1];1]
            new_point=new_point/new_point[4,1]
            if min>norm(original_point-new_point)
                min=norm(original_point-new_point)
            end
        end

        error_sum = error_sum+min
    end

   return log(10,(error_sum/s)/ε)
end

function Rec_L1_2_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspondences using (L1).2
    #
    
    errSum = []; 

    for i in 1:iter errSum=vcat(errSum,Rec_L1_2(m, s, ε)) end

    return errSum 
end

function Rec_L1_2_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average time of reconstructing s point correspondence using (L1).2
    #
    
    time = []; 

    for i in 1:iter time=vcat(time,@elapsed Rec_L1_2(m, s, ε)) end

    return time 
end


##

## (L1).3 ----------------------------------------------------------

##


∧(M) = M[:,1] × M[:,2];

function Rec_L1_3(m, s, ε)

    ## Input: m views, s point correspondences, ε noise level
    #
    ## Output: The average log relative error of reconstructing s point correspondences using (L1).3
    #

    cams = [randn(3, 4) for _ in 1:m]
    cams_vec = vec(cams[1]) 
    for t in 2:m cams_vec=vcat(cams_vec,vec(cams[t])) end

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0]; # Random line

    k = sort([randn(1) for i in 1:s]);
    PP = [line*[k[j];1] for j in 1:s]; # Points in R^3
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s]; # Images 
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s]; # Normalised images

    ll = [(∧(cams[i]*line)) for i in 1:m];
    ll = [ll[i]/ll[i][3] for i in 1:m];
    
    ## adding noise to points
    ee=rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    # Reconstruct the first point 
    qq_vec = qq[Int(floor((s+1)/2))][1][1:2] # We choose the middle point out of s many, since it might reduce the error
    for t in 2:m qq_vec = vcat(qq_vec,qq[Int(floor((s+1)/2))][t][1:2]) end 
    
    R = real_solutions(HomotopyContinuation.solve(D_PMV, start_PMV , start_parameters = PMV_q₀,  target_parameters = [cams_vec;qq_vec], show_progress=false))  

    
    original_point = PP[Int(floor((s+1)/2))]/PP[Int(floor((s+1)/2))][4,1]
    X_best = [R[1][1:3];1]
    min = norm(original_point - X_best)
    for k in 2:length(R)
        new_point = [R[k][1:3];1]
        if min > norm(original_point - new_point)
            min = norm(original_point - new_point)
            X_best = new_point
        end
    end
    
    error_sum = min

    X_plane = nullspace(transpose([X_best [0; 0; 0; 0]])) # ON-basis for a plane isomorphic to \Lambda(X_best)

    ee=rand(2,m)
    uu = [[(ll[i][1:2] + ε*ee[:,i]/norm(ee[:,i]) );1] for i in 1:m];
    
    AA_21 = [nullspace(transpose(cams[i]*[X_best [0; 0; 0; 0]])) for i in 1:m]
    uu_A = [transpose(AA_21[i])*[uu[i][1:2];1] for i in 1:m]
    uu_A = [uu_A[i]/uu_A[i][2,1] for i in 1:m]

    uu_A_vec = uu_A[1][1]
    for t in 2:m uu_A_vec=vcat(uu_A_vec,uu_A[t][1]) end
    
    cross_X = [[0 -(cams[i]*X_best)[3] (cams[i]*X_best)[2];
                (cams[i]*X_best)[3] 0 -(cams[i]*X_best)[1];
                -(cams[i]*X_best)[2] (cams[i]*X_best)[1] 0] for i in 1:m]
    cams_aug = [transpose(AA_21[i])*cross_X[i]*cams[i]*X_plane for i in 1:m]
    cams_aug_vec = vec(cams_aug[1])
    for t in 2:m cams_aug_vec=vcat(cams_aug_vec,vec(cams_aug[t])) end 

    R = real_solutions(HomotopyContinuation.solve(D_M21, start_M21 , start_parameters = M21_q₀,  target_parameters = [cams_aug_vec;uu_A_vec],show_progress=false))  

    min = 10
    L_best = 0 
    for k in 1:length(R)
        new_point = X_plane*[R[k][1:2];1]
        new_point = new_point/new_point[4,1]
        new_line = hcat(X_best,new_point)
        new_line = new_line*(inv(new_line[3:4,:]))
        original_line = line*(inv(line[3:4,:]))
        if min > norm(original_line - new_line)
            min = norm(original_line - new_line)
            L_best = new_line
        end
    end

    AA = [transpose(Matrix(qr(cams[i]*L_best).Q)) for i in 1:m]  # Projection matrices X such that XX^T=I  
    qq_A = [[AA[i]*qq[j][i] for i in 1:m] for j in 1:s] 
    qq_A = [[qq_A[j][i]/qq_A[j][i][2] for i in 1:m] for j in 1:s]

    cams_aug = [AA[i]*cams[i]*L_best for i in 1:m]
    cams_aug_vec = vcat(vec(cams_aug[1]))
    for t in 2:m cams_aug_vec=vcat(cams_aug_vec,vec(cams_aug[t])) end 
 
    for j in 1:s # Reconstruct all the point correspondences in the anchored PMV
        if j != Int(floor((s+1)/2))

            qq_A_vec = qq_A[j][1][1]
            for t in 2:m qq_A_vec=vcat(qq_A_vec, qq_A[j][t][1]) end

            R = real_solutions(HomotopyContinuation.solve(D_M11, start_M11 , start_parameters = M11_q₀, 
                                                          target_parameters = [cams_aug_vec;qq_A_vec],show_progress=false))
            
            min=10
            original_point=PP[j]/PP[j][4,1]
            for k in 1:length(R)
                new_point=L_best*[R[k][1];1]
                new_point=new_point/new_point[4,1]
                if min>norm(original_point-new_point)
                    min=norm(original_point-new_point)
                end
            end

            error_sum = error_sum+min
        end
    end

    return log(10,(error_sum/s)/ε)
end


function Rec_L1_3_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspondences using (L1).3
    #

    errSum = [];
    
    for i in 1:iter errSum=vcat(errSum,Rec_L1_3(m,s,ε)) end

    return errSum  
end

function Rec_L1_3_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average time of reconstructing s point correspondences using (L1).3
    #

    time = [];
    
    for i in 1:iter time=vcat(time,@elapsed Rec_L1_3(m,s,ε)) end

    return time  
end


##

## (L1).4 ----------------------------------------------------------

##

∧(M) = M[:,1] × M[:,2];

function Rec_L1_4(m, s, ε)

    ## Input: m views, s point correspondences, ε noise level
    #
    ## Output: The average log relative error of reconstructing s point correspondences using (L1).4
    #

    cams = [randn(3, 4) for _ in 1:m]
    cams_vec = vec(cams[1]) 
    for t in 2:m cams_vec=vcat(cams_vec,vec(cams[t])) end

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];1 0;0 1] # Random line

    k = randn(1,s);
    PP = [line*[k[j];1] for j in 1:s]; # Points in R^3
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s]; # Images 
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s]; # Normalised images

    ll = [(∧(cams[i]*line)) for i in 1:m];
    ll = [ll[i]/ll[i][3] for i in 1:m];
    
    ## adding noise to points
    ee=rand(s,2,m)
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    ee = rand(2,m)
    uu = [[(ll[i][1:2] + ε*ee[:,i]/norm(ee[:,i]) );1] for i in 1:m];
    
    uu_vec = uu[1][1:2]
    for i in 2:m uu_vec=vcat(uu_vec,uu[i][1:2]) end

    R = real_solutions(HomotopyContinuation.solve(D_LMV, start_LMV , start_parameters = LMV_q₀,  target_parameters = [cams_vec;uu_vec],show_progress=false))
            
    min = 10
    L_best = 0 # Initialize the best reconstructed line
    for k in 1:length(R)
        new_line = [R[k][1] R[k][3];R[k][2] R[k][4];1 0;0 1]
        original_line = line 
        if min > norm(original_line - new_line) # Approximates the distance of lines
            min = norm(original_line - new_line)
            L_best = new_line
        end
    end

    AA = [transpose(Matrix(qr(cams[i]*L_best).Q)) for i in 1:m]  # Projection matrices X such that XX^T=I  
    qq_A = [[AA[i]*qq[j][i] for i in 1:m] for j in 1:s] 
    qq_A = [[qq_A[j][i]/qq_A[j][i][2] for i in 1:m] for j in 1:s]

    cams_aug = [AA[i]*cams[i]*L_best for i in 1:m]
    cams_aug_vec = vcat(vec(cams_aug[1]))
    for t in 2:m cams_aug_vec=vcat(cams_aug_vec,vec(cams_aug[t])) end 
    
    error_sum = 0

    for j in 1:s #

        qq_A_vec = qq_A[j][1][1]
        for t in 2:m qq_A_vec=vcat(qq_A_vec, qq_A[j][t][1]) end

        R = real_solutions(HomotopyContinuation.solve(D_M11, start_M11 , start_parameters = M11_q₀,  target_parameters = [cams_aug_vec;qq_A_vec],show_progress=false))
        
        min=10
        original_point=PP[j]/PP[j][4,1]
        for k in 1:length(R)
            new_point=L_best*[R[k][1];1]
            new_point=new_point/new_point[4,1]
            if min>norm(original_point-new_point)
                min=norm(original_point-new_point)
            end
        end

        error_sum = error_sum+min

    end
    return log(10,(error_sum/s)/ε)
end


function Rec_L1_4_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspondences using (L1).4
    #

    errSum = [];
    
    for i in 1:iter errSum=vcat(errSum,Rec_L1_4(m,s,ε)) end

    return errSum  
end

function Rec_L1_4_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average time of reconstructing s point correspondences using (L1).4
    #

    time = [];
    
    for i in 1:iter time=vcat(time,@elapsed Rec_L1_4(m,s,ε)) end

    return time  
end

##

## Comparison between standard and non-standard methods

##


##

## (L1).1 std ----------------------------------------------------------

##



function Rec_L1_1_std(m, s, ε)
    
    ## Input: m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspodences using (L1).1 standard approach
    #

    cams = [randn(3, 4) for _ in 1:m];
    cams_vec=vec(cams[1])
    for t in 2:m cams_vec=vcat(cams_vec,vec(cams[t])) end

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0]; # Random line

    k = randn(1, s);
    PP = [line*[k[j];1] for j in 1:s]; #points
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s];
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s];
    
    ll = [(∧(cams[i]*line)) for i in 1:m];
    ll = [ll[i]/ll[i][3] for i in 1:m];

    ## adding noise to points
    ee=rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    ee=rand(2,m)
    uu = [[(ll[i][1:2] + ε*ee[:,i]/norm(ee[:,i]) );1] for i in 1:m];
    
    bp_line = nullspace(transpose([transpose(cams[1])*uu[1] transpose(cams[2])*uu[2]]));
    LL_vec = vec(bp_line)

    error_sum = 0
    
    for j in 1:s # Reconstruct all the point correspondences in the anchored PMV

        qq_vec = qq[j][1][1:2]
        for t in 2:m qq_vec=vcat(qq_vec,qq[j][t][1:2]) end
        R = real_solutions(HomotopyContinuation.solve(D_aPMV, start_aPMV , start_parameters = aPMV_q₀,  target_parameters = [cams_vec;qq_vec;LL_vec],show_progress=false))
        
        min=10
        original_point=PP[j]/PP[j][4,1]
        for k in 1:length(R)
            new_point = bp_line*[R[k][1];1]
            new_point = new_point/new_point[4,1]
            if min>norm(original_point-new_point)
                min=norm(original_point-new_point)
            end
        end

        error_sum = error_sum+min
    end

   return log(10,(error_sum/s)/ε)
end

function Rec_L1_1_std_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing using (L1).1 standard approach
    #

    errSum = [];
    
    for i in 1:iter errSum=vcat(errSum,Rec_L1_1_std(m,s,ε)) end

    return errSum  
end

function Rec_L1_1_std_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average time of reconstructing using (L1).1 standard approach
    #

    time = [];
    
    for i in 1:iter time=vcat(time,@elapsed Rec_L1_1_std(m,s,ε)) end

    return time  
end


##

## (L1).3 std ----------------------------------------------------------

##


function Rec_L1_3_std(m, s, ε)

    ## Input: m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing s point correspodences using (L1).1 standard approach
    #

    cams = [randn(3, 4) for _ in 1:m];
    cams_vec = vec(cams[1])
    for t in 2:m cams_vec = vcat(cams_vec,vec(cams[t])) end

    v = randn(2, 2); 
    line = [v[1,1] v[1,2];v[2,1] v[2,2];0 1;1 0]; # Random line

    
    k = sort([randn(1) for i in 1:s]);
    PP = [line*[k[j];1] for j in 1:s]; # Points in R^3
    pp = [[cams[i]*PP[j] for i in 1:m] for j in 1:s]; # Images 
    pp = [[pp[j][i]/pp[j][i][3] for i in 1:m] for j in 1:s]; # Normalised images

    ll = [(∧(cams[i]*line)) for i in 1:m];
    ll = [ll[i]/ll[i][3] for i in 1:m];
    
    ## adding noise to points
    ee=rand(s,2,m) 
    qq = [[[(pp[j][i][1:2] + ε*ee[j,:,i]/norm(ee[j,:,i]) );1] for i in 1:m] for j in 1:s];

    # Reconstruct the first point 
    qq_vec = qq[Int(floor((s+1)/2))][1][1:2] # We choose the middle point out of s many, since it might reduce the error
    for t in 2:m qq_vec = vcat(qq_vec,qq[Int(floor((s+1)/2))][t][1:2]) end 
    
    R = real_solutions(HomotopyContinuation.solve(D_PMV, start_PMV , start_parameters = PMV_q₀,  target_parameters = [cams_vec;qq_vec], show_progress=false))  

    original_point = PP[Int(floor((s+1)/2))]/PP[Int(floor((s+1)/2))][4,1]
    X_best = [R[1][1:3];1]
    min = norm(original_point - X_best)
    for k in 2:length(R)
        new_point = [R[k][1:3];1]
        if min > norm(original_point - new_point)
            min = norm(original_point - new_point)
            X_best = new_point
        end
    end
    
    error_sum = min

    ee=rand(2,m)
    uu = [[(ll[i][1:2] + ε*ee[:,i]/norm(ee[:,i]) );1] for i in 1:m];
    uu_vec = uu[1][1:2]
    for t in 2:m uu_vec = vcat(uu_vec,uu[t][1:2]) end
    #
    R = real_solutions(HomotopyContinuation.solve(D_aLMV, start_aLMV , start_parameters = aLMV_q₀,  target_parameters = [X_best[1:3];cams_vec;uu_vec],show_progress=false))  

    min=10
    L_best=0 
    for k in 1:length(R)
        new_line = hcat(X_best,[R[k][1:2];1;0])
        new_line = hcat(new_line[1:4,1]-new_line[3,1]*new_line[1:4,2], new_line[1:4,2]) 
        if min>norm(line-new_line)
            min=norm(line-new_line)
            L_best=new_line
        end
    end

    LL_vec = vec(L_best)
    
    for j in 1:s # Reconstruct all the point correspondences in the anchored PMV
        if j != Int(floor((s+1)/2))

            qq_vec = qq[j][1][1:2]
            for t in 2:m qq_vec = vcat(qq_vec,qq[j][t][1:2]) end
            R = real_solutions(HomotopyContinuation.solve(D_aPMV, start_aPMV , start_parameters = aPMV_q₀,  target_parameters = [cams_vec;qq_vec;LL_vec],show_progress=false))
            
            min=10
            original_point = PP[j]/PP[j][4,1]
            for k in 1:length(R)
                new_point = L_best*[R[k][1];1]
                new_point = new_point/new_point[4,1]
                if min > norm(original_point-new_point)
                    min = norm(original_point-new_point)
                end
            end

            error_sum = error_sum+min
        end
    end

   return log(10,(error_sum/s)/ε)
end


function Rec_L1_3_std_Error(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing using (L1).3 standard approach
    #

    errSum = [];
    
    for i in 1:iter errSum=vcat(errSum,Rec_L1_3_std(m,s,ε)) end

    return errSum  
end

function Rec_L1_3_std_Time(iter, m, s, ε)
    
    ## Input: iter iterations, m cameras, s point correspondences, ε noise level
    #
    ## Output: The average error of reconstructing using (L1).3 standard approach
    #

    time = [];
    
    for i in 1:iter time=vcat(time,@elapsed Rec_L1_3_std(m,s,ε)) end

    return time 
end

