using HomotopyContinuation;
using LinearAlgebra;

## Noise Level
ε = 1e-12;
#

m = 4; # Number of views
s = 5; # Number of points correspondences incident to our line correspondence

##

## Denoising in the PMV ----------------------------------------------------------------

##

@var X[1:3]; # Variable: Point in R^3 
@var cams[1:m,1:3,1:4] # Parameter: Cameras
@var qq[1:m,1:2]; # Parameter: Noisy images
@var t[1:m]; # Variable: Projective scaling
@var λ1[1:m]; # Variable: Lagrange multipliers

pp = [[(t[i]*cams[i,1:3,1:4]*[X;1])[1:2];1] for i in 1:m]; # X projected onto the different camera planes
dif = [pp[i] - [qq[i,1:2];1] for i in 1:m]; # Difference between projected points and noisy image points 
F = transpose(dif) * (dif) 
U = sum(λ1[i] * ((t[i]*cams[i,1:3,1:4]*[X;1])[3]-1) for i in 1:m) # Lagrange multiplier condition
L = F - U # Lagrange function
∇L = differentiate(L, [X;t;λ1]) 

cams_vec = vec(cams[1,:,:]) #vcat(cams[1,1:3,1],cams[1,1:3,2],cams[1,1:3,3],cams[1,1:3,4]) # Stacking camera variables  
qq_vec = vec(qq[1,:]) #qq[1,1:2] # Stacking image point variables
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))#cams[i,1:3,1],cams[i,1:3,2],cams[i,1:3,3],cams[i,1:3,4])
             global qq_vec=vcat(qq_vec,vec(qq[i,:]))#qq[i,1:2])  
end

D_PMV = System(∇L; variables = [X;t;λ1], parameters = [cams_vec;qq_vec]) # Polynomial system


## Initial Step
PMV_q₀ = randn(ComplexF64, 12*m+2*m) # Random starting parameters
PMV_S₀ = solve(D_PMV; target_parameters = PMV_q₀) # Starting system
start_PMV = solutions(PMV_S₀) # Starting solution

##

## Denoising in M^{1,1} --------------------------------------------------

##


@var X; # Variable: Point in R^1 
@var cams[1:m,1:2,1:2]; # Parameter: Cameras 
@var qq[1:m,1]; # Parameter: Noisy images
@var t[1:m]; # Variable: Projective scaling
@var λ1[1:m]; # Variable: Lagrange multipliers

pp = [[(t[i]*cams[i,1:2,1:2]*[X;1])[1];1] for i in 1:m]; # Projected points
dif = [pp[i] - [qq[i,1];1] for i in 1:m];
F = transpose(dif) * (dif)
U = sum(λ1[i] * ((t[i]*cams[i,1:2,1:2]*[X;1])[2]-1) for i in 1:m)
L = F - U
∇L = differentiate(L, [X;t;λ1])

cams_vec = vec(cams[1,:,:])
qq_vec = qq[1,1]
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))#cams[i,1:2,1],cams[i,1:2,2])
             global qq_vec=vcat(qq_vec,qq[i,1])  
end

D_M11 = System(∇L; variables = [X;t;λ1], parameters = [cams_vec;qq_vec])


## Initial Step
M11_q₀ = randn(ComplexF64, 4*m+m)
M11_S₀ = solve(D_M11; target_parameters = M11_q₀)
start_M11 = solutions(M11_S₀)


##

## Denoising in M^{2,1} ---------------------------------------------

##


@var X[1:2]; # Variable: Point in R^3 
@var cams[1:m,1:2,1:3] # Parameter: Cameras
@var qq[1:m,1]; # Parameter: Noisy images
@var t[1:m]; # Variable: Projective scaling
@var λ1[1:m]; # Variable: Lagrange multipliers

pp = [[(t[i]*cams[i,1:2,1:3]*[X;1])[1];1] for i in 1:m]; # Projected points
dif = [pp[i] - [qq[i,1];1] for i in 1:m];
F = transpose(dif) * (dif)
U = sum(λ1[i] * ((t[i]*cams[i,1:2,1:3]*[X;1])[2]-1) for i in 1:m)
L = F - U
∇L = differentiate(L, [X;t;λ1])

cams_vec = vec(cams[1,:,:])
qq_vec = qq[1,1]
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))
             global qq_vec=vcat(qq_vec,qq[i,1])  
end

D_M21 = System(∇L; variables = [X;t;λ1], parameters = [cams_vec;qq_vec])


## Initial Step
M21_q₀ = randn(ComplexF64, 6*m+m)
M21_S₀ = solve(D_M21; target_parameters = M21_q₀)
start_M21 = solutions(M21_S₀)


##

## Denoising in LMV ---------------------------------------------

##



∧(M) = M[:,1] × M[:,2];

@var v[1:2,1:2]
@var cams[1:m,1:3,1:4]
@var uu[1:m,1:2]
@var p[1:m]
@var λ1[1:m]

L = [v[1,1] v[1,2];v[2,1] v[2,2];1 0;0 1]
l = [[(p[i].*∧(cams[i,1:3,1:4]*L))[1:2];1] for i in 1:m]
g = [l[i] - [uu[i,1:2];1] for i in 1:m]
F = transpose(g) * (g)

U = sum(λ1[i] * ((p[i].*∧(cams[i,1:3,1:4]*L))[3]-1) for i in 1:m)
L = F - U
vars = [vec(v);p;λ1] 
∇L = differentiate(L, vars)


cams_vec = vec(cams[1,:,:]) 
uu_vec = vec(uu[1,:]) #qq[1,1:2] # Stacking image point variables
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))
             global uu_vec=vcat(uu_vec,vec(uu[i,:])) 
end


D_LMV = System(∇L; variables = vars, parameters = [cams_vec;uu_vec])

LMV_q₀ = randn(ComplexF64, 12*m+2*m) # Random starting parameters
LMV_S₀ = solve(D_LMV; target_parameters = LMV_q₀, show_progress=true) # Starting system
start_LMV = solutions(LMV_S₀) # Starting solution



##

## Denoising in the anchored PMV ------------------------------------------------------------------------------

##


@var v; # Variable: Point on the starting line
@var cams[1:m,1:3,1:4]; # Parameter: Cameras
@var qq[1:m,1:2]; # Parameter: Noisy images 
@var LL[1:4,1:2]; # Parameter: Line 
@var t[1:m]; # Variable: Projective scaling
@var λ1[1:m]; # Variable: Lagrange multipliers

PP=LL*[v;1] # The point determined by v
pp = [[(t[i]*cams[i,1:3,1:4]*PP)[1:2];1] for i in 1:m]; # Projected points
dif = [pp[i] - [qq[i,1:2];1] for i in 1:m];
F = transpose(dif) * (dif)
U = sum(λ1[i] * ((t[i]*cams[i,1:3,1:4]*PP)[3]-1) for i in 1:m)
L = F - U
∇L = differentiate(L, [v;t;λ1])

LL_vec = vec(LL) # vcat(LL[1:4,1],LL[1:4,2]) # Stacking line variables
cams_vec = vec(cams[1,:,:]) # Stacking camera variables 
qq_vec = qq[1,1:2] # Stacking noisy image variables
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))
             global qq_vec=vcat(qq_vec,qq[i,1:2])  
end

D_aPMV = System(∇L; variables = [v;t;λ1], parameters = [cams_vec;qq_vec;LL_vec])

## Initial Step  ## Improve by not choosing randomg Anc_q0
aPMV_q₀ = randn(ComplexF64, 12*m+2*m+8)
aPMV_S₀ = solve(D_aPMV; target_parameters = aPMV_q₀) 
start_aPMV = solutions(aPMV_S₀)



##

## Denoising in the anchored LMV ------------------------------------------------------------------------------

##

∧(M) = M[:,1] × M[:,2]; # Describes the joint image map for lines

@var X[1:3]; # Parameter: The fixed point X determining the anchored LMV
@var v[1:2]; # Variable: Other point of the line
@var cams[1:m,1:3,1:4]; # Parameter: Cameras
@var uu[1:m,1:2]; # Parameter: noisy lines
@var t[1:m]; # Variable: Projective scaling
@var λ1[1:m]; # Variable: Lagrange multipliers

LL=hcat([X;1],[v;1;0]) # The line determine by X and v
ll = [[(t[i]*∧(cams[i,1:3,1:4]*LL))[1:2];1] for i in 1:m]; # Projected points
dif = [ll[i] - [uu[i,1:2];1] for i in 1:m];
F = transpose(dif) * (dif)
U = sum(λ1[i] * ((t[i]*∧(cams[i,1:3,1:4]*LL))[3]-1) for i in 1:m)
L = F - U
∇L = differentiate(L, [v;t;λ1])

cams_vec = vec(cams[1,:,:])
uu_vec = uu[1,1:2]
for i in 2:m global cams_vec=vcat(cams_vec,vec(cams[i,:,:]))
             global uu_vec=vcat(uu_vec,uu[i,1:2])  
end

D_aLMV = System(∇L; variables = [v;t;λ1], parameters = [X;cams_vec;uu_vec])

aLMV_q₀ = randn(ComplexF64, 3+12*m+2*m)
aLMV_S₀ = solve(D_aLMV; target_parameters = aLMV_q₀)
start_aLMV = solutions(aLMV_S₀)
