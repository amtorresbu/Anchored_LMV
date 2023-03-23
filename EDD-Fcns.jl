using HomotopyContinuation;
using LinearAlgebra;

# Functions numerical evidence of EDDs of multiview varieties from the main body of the paper


function EDD_aPMV(m) 

    ## Input: m number of cameras
    #
    ## Output: Number of critical points for the ED problem given m random cameras and random image points
    #

    cams = [randn(3, 4) for _ in 1:m]

    @var v
    @var p[1:m] 
    @var λ1[1:m]

    PP = [v; 0; 0; 1]
    pp = [[(p[i]*cams[i]*PP)[1:2];1] for i in 1:m]
    qq = randn(2, m)
    dif = [pp[i] - [qq[:,i];1] for i in 1:m]
    F = transpose(dif) * (dif)
    U = sum(λ1[i] * ((p[i]*cams[i]*PP)[3]-1) for i in 1:m)
    L = F - U
    ∇L = differentiate(L, [v;p;λ1] )
    D = System(∇L; variables = [v;p;λ1] )
    R = solutions(HomotopyContinuation.solve(D, show_progress=true))
    nonSing = filter(r->all(abs.(r).>1e-8), R)
    
    O = certify(D,nonSing)
    o = certificates(O)
    lll = map(o) do oᵢ
        if !is_certified(oᵢ)
            return false
        end
        𝕀 = certified_solution_interval(oᵢ)
        for i in 1:m
            if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                return false
            end
        end
        return true
    end

    EDD = length(unique_points(nonSing)), count(lll)

    return println("Number of critical points for $m cameras: $EDD, Expected: ", Int(3*m-2))

end
        

function EDD_aLMV(m) 

    ## Input: m number of cameras
    #
    ## Output: Number of critical points for the ED problem given m random cameras and random image lines
    #

    cams = [randn(3, 4) for _ in 1:m]
    ∧(M) = M[:,1] × M[:,2]

    @var v[1:2]
    @var p[1:m]
    @var λ1[1:m]

    LL = [1 0; v[1] 0; v[2] 0; 0 1]
    ll = [[(p[i].*∧(cams[i]*LL))[1:2];1] for i in 1:m]
    uu = randn(2, m)
    dif = [ll[i] - [uu[:,i];1] for i in 1:m]
    F = transpose(dif) * (dif)    
    U = sum(λ1[i] * ((p[i].*∧(cams[i]*LL))[3]-1) for i in 1:m)
    L = F - U 
    ∇L = differentiate(L, [vec(v);p;λ1])
    D = System(∇L; variables = [vec(v);p;λ1])
    R = solutions(HomotopyContinuation.solve(D, show_progress=true))
    nonSing = filter(r->all(abs.(r).>1e-8), R)
    
    O = certify(D,nonSing)
        o = certificates(O)
        lll = map(o) do oᵢ
            if !is_certified(oᵢ)
                return false
            end
            𝕀 = certified_solution_interval(oᵢ)
            for i in 1:m
                if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                    return false
                end
            end
            return true
        end
    
    EDD = length(unique_points(nonSing)), count(lll)

    return println("Number of critical points for $m cameras: $EDD, Expected: ", Int((9/2)*m^2-(19/2)*m+3))
    
end

        
function EDD_LMV(m) 

    ## Input: m number of cameras
    # 
    ## Output: Number of critical points for the ED problem given m random cameras and random image lines
    #

    cams = [randn(3, 4) for _ in 1:m]
    ∧(M) = M[:,1] × M[:,2]

    @var v[1:2,1:2]
    @var p[1:m]
    @var λ1[1:m]

    LL = [1 0; 0 1; v]
    ll = [[(p[i].*∧(cams[i]*LL))[1:2];1] for i in 1:m]
    uu = randn(2, m)
    dif = [ll[i] - [uu[:,i];1] for i in 1:m]
    F = transpose(dif) * (dif)    
    U = sum(λ1[i] * ((p[i].*∧(cams[i]*LL))[3]-1) for i in 1:m)
    L = F - U 
    ∇L = differentiate(L, [vec(v);p;λ1])
    D = System(∇L; variables = [vec(v);p;λ1])
    R = solutions(HomotopyContinuation.solve(D, show_progress=true))
    nonSing = filter(r->all(abs.(r).>1e-8), R)
    
    O = certify(D,nonSing)
        o = certificates(O)
        lll = map(o) do oᵢ
            if !is_certified(oᵢ)
                return false
            end
            𝕀 = certified_solution_interval(oᵢ)
            for i in 1:m
                if Base.in(0, HomotopyContinuation.IComplexF64(𝕀[i]))
                    return false
                end
            end
            return true
        end

    EDD = length(unique_points(nonSing)), count(lll)

    return println("Number of critical points for $m cameras: $EDD")
end
