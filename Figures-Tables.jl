
## ----------------------------------------------------------------------

## First Plot

## Comparison (L1).1 and (L1).3 using standard and non-standard approaches

iter = 1000

## Load Start-Systems.jl with m = 2 and Reconstruction-Fcns.jl

# Note that if m = 2, then (L1).0 ~ (L1).2 and (L1).1 ~ (L1).4

error_L1_0 = Rec_L1_0_Error(iter,m,s,ε)
error_L1_1 = Rec_L1_1_Error(iter,m,s,ε)
error_L1_3 = Rec_L1_3_Error(iter,m,s,ε) # Iterations might fail a few times, try running the function for iter = 200 and concatening 5 times
error_L1_1_std = Rec_L1_1_std_Error(iter,m,s,ε) 
error_L1_3_std = Rec_L1_3_std_Error(iter,m,s,ε) # Iterations might fail a few times, try running the function for iter = 200 and concatening 5 times


time_L1_0 = Rec_L1_0_Time(iter,m,s,ε)
time_L1_1 = Rec_L1_1_Time(iter,m,s,ε)
time_L1_3 = Rec_L1_3_Time(iter,m,s,ε) # Iterations might fail a few times, try running the function for iter = 200 and concatening 5 times
time_L1_1_std = Rec_L1_1_std_Time(iter,m,s,ε)
time_L1_3_std = Rec_L1_3_std_Time(iter,m,s,ε) # Iterations might fail a few times, try running the function for iter = 200 and concatening 5 times

using Statistics

median(error_L1_0)
median(error_L1_1)
median(error_L1_3)
median(error_L1_1_std)
median(error_L1_3_std)


mean(error_L1_0)
mean(error_L1_1)
mean(error_L1_3)
mean(error_L1_1_std)
mean(error_L1_3_std)


std(error_L1_0)
std(error_L1_1)
std(error_L1_3)
std(error_L1_1_std)
std(error_L1_3_std)


median(time_L1_0[2:end])
median(time_L1_1[2:end])
median(time_L1_3[2:end])
median(time_L1_1_std[2:end])
median(time_L1_3_std[2:end])


mean(time_L1_0[2:end])
mean(time_L1_1[2:end])
mean(time_L1_3[2:end])
mean(time_L1_1_std[2:end])
mean(time_L1_3_std[2:end])


std(time_L1_0[2:end])
std(time_L1_1[2:end])
std(time_L1_3[2:end])
std(time_L1_1_std[2:end])
std(time_L1_3_std[2:end])


using Plots;


##

## Plot: Comparison of all (L1).0-4 methods ------------------------

##


# Accuracy:

fig_error = histogram(error_L1_0,
            bins = -2:0.1:6,
            alpha = 0.6,
            label = "(L1).0",
            title = "Accuracy: $m cameras, $s point correspondences")

histogram!(error_L1_1_std,
            bins = -2:0.1:6,
            color = "grey",
            alpha = 0.6,
            label = "(L1).1 std")

histogram!(error_L1_3_std,
            bins = -2:0.1:6,
            color = "brown",
            alpha = 0.6,
            label = "(L1).3 std")
histogram!(error_L1_1,
            bins = -2:0.1:6,
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(error_L1_3,
            bins = -2:0.1:6,
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")

xlabel!("Relative error")
ylabel!("Frequency")
savefig(fig_error, "Hist-error-1000-std-nonstd-$m-cam.png")

# Speed:

fig_time = histogram(time_L1_0,
            bins = 0:0.0025:0.3,#0:0.008:0.8, # 0:0.0025:0.3,
            alpha = 0.6,
            label = "(L1).0",
            title = "Speed: $m cameras, $s point correspondences")
histogram!(time_L1_1_std,
            bins = 0:0.0025:0.3,#0:0.008:0.8
            color = "grey",
            alpha = 0.6,
            label = "(L1).1 std")
histogram!(time_L1_3_std,
            bins = 0:0.0025:0.3,#0:0.008:0.8
            color = "brown",
            alpha = 0.6,
            label = "(L1).3 std")
histogram!(time_L1_1,
            bins = 0:0.0025:0.3,#0:0.008:0.8
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(time_L1_3,
            bins = 0:0.0025:0.3,#0:0.008:0.8
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")
xlabel!("Seconds")
ylabel!("Frequency")
savefig(fig_time, "Hist-time-1000-std-nonstd-$m-cam.png")


## ----------------------------------------------------------------------

## Second Plot

## Comparison of (L1).0-4 

## Load Start-Systems.jl with m = 3 and Reconstruction-Fcns.jl

iter = 1000

error_L1_0 = Rec_L1_0_Error(iter,m,s,ε)
error_L1_1 = Rec_L1_1_Error(iter,m,s,ε)
error_L1_2 = Rec_L1_2_Error(iter,m,s,ε)
error_L1_3 = Rec_L1_3_Error(iter,m,s,ε) 
error_L1_4 = Rec_L1_4_Error(iter,m,s,ε)


time_L1_0 = Rec_L1_0_Time(iter,m,s,ε)
time_L1_1 = Rec_L1_1_Time(iter,m,s,ε)
time_L1_2 = Rec_L1_2_Time(iter,m,s,ε)
time_L1_3 = Rec_L1_3_Time(iter,m,s,ε) 
time_L1_4 = Rec_L1_4_Time(iter,m,s,ε)


using Statistics

median(error_L1_0)
median(error_L1_1)
median(error_L1_2)
median(error_L1_3)
median(error_L1_4)

mean(error_L1_0)
mean(error_L1_1)
mean(error_L1_2)
mean(error_L1_3)
mean(error_L1_4)

std(error_L1_0)
std(error_L1_1)
std(error_L1_2)
std(error_L1_3)
std(error_L1_4)

median(time_L1_0)
median(time_L1_1)
median(time_L1_2)
median(time_L1_3)
median(time_L1_4)

mean(time_L1_0[2:1000])
mean(time_L1_1[2:1000])
mean(time_L1_2[2:1000])
mean(time_L1_3[2:1000])
mean(time_L1_4[2:1000])

std(time_L1_0[2:1000])
std(time_L1_1[2:1000])
std(time_L1_2[2:1000])
std(time_L1_3[2:1000])
std(time_L1_4[2:1000])

using Plots;


##

## Plot: Comparison of all (L1).0-4 methods ------------------------

##


# Accuracy:

fig_error = histogram(error_L1_0,
            bins = -2:0.2:6,
            alpha = 0.6,
            label = "(L1).0",
            title = "Accuracy: $m cameras, $s point correspondences")
histogram!(error_L1_1,
            bins = -2:0.2:6,
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(error_L1_2,
            bins = -2:0.2:6,
            color = "orange",
            alpha = 0.6,
            label = "(L1).2")
histogram!(error_L1_3,
            bins = -2:0.2:6,
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")
histogram!(error_L1_4,
            bins = -2:0.2:6,
            color = "cyan",
            alpha = 0.6,
            label = "(L1).4")
xlabel!("Relative error")
ylabel!("Frequency")
savefig(fig_error, "Hist-error-1000-01234-$m-cam.png")

# Speed:

fig_time = histogram(time_L1_0,
            bins = 0:0.01875:1.5, 
            alpha = 0.6,
            label = "(L1).0",
            title = "Speed: $m cameras, $s point correspondences")

histogram!(time_L1_1,
            bins =  0:0.01875:1.5,
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(time_L1_2,
            bins =  0:0.01875:1.5,
            color = "orange",
            alpha = 0.6,
            label = "(L1).2")
histogram!(time_L1_3,
            bins =  0:0.01875:1.5,
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")
histogram!(time_L1_4,
            bins = 0:0.01875:1.5,
            color = "cyan",
            alpha = 0.6,
            label = "(L1).4")
xlabel!("Seconds")
ylabel!("Frequency")
savefig(fig_time, "Hist-time-1000-01234-$m-cam.png")


## ----------------------------------------------------------------------

## Third Plot

## Comparison of (L1).0-4 

## Load Start-Systems.jl with m = 4 and Reconstruction-Fcns.jl

iter = 100

error_L1_0 = Rec_L1_0_Error(iter,m,s,ε)
error_L1_1 = Rec_L1_1_Error(iter,m,s,ε)
error_L1_2 = Rec_L1_2_Error(iter,m,s,ε)
error_L1_3 = Rec_L1_3_Error(iter,m,s,ε) 
error_L1_4 = Rec_L1_4_Error(iter,m,s,ε)


time_L1_0 = Rec_L1_0_Time(iter,m,s,ε)
time_L1_1 = Rec_L1_1_Time(iter,m,s,ε)
time_L1_2 = Rec_L1_2_Time(iter,m,s,ε)
time_L1_3 = Rec_L1_3_Time(iter,m,s,ε) 
time_L1_4 = Rec_L1_4_Time(iter,m,s,ε)


using Statistics

median(error_L1_0)
median(error_L1_1)
median(error_L1_2)
median(error_L1_3)
median(error_L1_4)

mean(error_L1_0)
mean(error_L1_1)
mean(error_L1_2)
mean(error_L1_3)
mean(error_L1_4)

std(error_L1_0)
std(error_L1_1)
std(error_L1_2)
std(error_L1_3)
std(error_L1_4)

median(time_L1_0)
median(time_L1_1)
median(time_L1_2)
median(time_L1_3)
median(time_L1_4)

mean(time_L1_0[2:100])
mean(time_L1_1[2:100])
mean(time_L1_2[2:100])
mean(time_L1_3[2:100])
mean(time_L1_4[2:100])

std(time_L1_0[2:100])
std(time_L1_1[2:100])
std(time_L1_2[2:100])
std(time_L1_3[2:100])
std(time_L1_4[2:100])

using Plots;

##

## Plot: Comparison of all (L1).0-4 methods ------------------------

##

# Accuracy:

fig_error = histogram(error_L1_0,
            bins = -2:0.2:6,
            alpha = 0.6,
            label = "(L1).0",
            title = "Accuracy: $m cameras, $s point correspondences")
histogram!(error_L1_1,
            bins = -2:0.2:6,
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(error_L1_2,
            bins = -2:0.2:6,
            color = "orange",
            alpha = 0.6,
            label = "(L1).2")
histogram!(error_L1_3,
            bins = -2:0.2:6,
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")
histogram!(error_L1_4,
            bins = -2:0.2:6,
            color = "cyan",
            alpha = 0.6,
            label = "(L1).4")
xlabel!("Relative error")
ylabel!("Frequency")
savefig(fig_error, "Hist-error-100-01234-$m-cam.png")

# Speed:

fig_time = histogram(time_L1_0,
            bins = 0:0.25:10,  
            alpha = 0.6,
            label = "(L1).0",
            title = "Speed: $m cameras, $s point correspondences")

histogram!(time_L1_1,
            bins =  0:0.25:10,
            color = "magenta",
            alpha = 0.6,
            label = "(L1).1")
histogram!(time_L1_2,
            bins =  0:0.25:10, 
            color = "orange",
            alpha = 0.6,
            label = "(L1).2")
histogram!(time_L1_3,
            bins =  0:0.25:10, 
            color = "yellow",
            alpha = 0.6,
            label = "(L1).3")
histogram!(time_L1_4,
            bins = 0:0.25:10,
            color = "cyan",
            alpha = 0.6,
            label = "(L1).4")
xlabel!("Seconds")
ylabel!("Frequency")
savefig(fig_time, "Hist-time-100-01234-$m-cam.png")
