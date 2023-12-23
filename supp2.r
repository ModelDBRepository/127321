#####################################################
#
# Simulation of the C.elegans egg-laying circuit
#
#####################################################
#
# Zhang Mi, William R. Schafer, R. Breitling
#
#####################################################


# Parameters 
timepoints = 12000
tau = 0.5						# size of simulation time step
N = 10; 						# replicates of stochastic simulation to be run

# Constants
Thalf=140  
p0=exp(-tau*log(2)/Thalf)			# half time of clearance is about 140 seconds.
lambda1 = 1/23*tau;
lambda2 = 1/1800*tau;
clusterN = 3                              # threshold indicating uv1 sensitivity to neurotransmitter and 
                                          # thus controlling the number of eggs per cluster

# Initial state transition probabilities
p1=1
p2=lambda1
p3=0
p4=lambda2
p5=0
p6=0
p7=1

eggnumber = CN = c();      			# output statistics: eggnumber = rate egg-laying events/second; CN, egg-laying events/cluster


for (k in 1:N) {

# Variable initialization
egg = eggn = 0
hsn = hsnn = 1
vc = vcn = 0
uv1 = uv1n = 0
count = 0
eggs = vcs = uv1s = counts = clusters = hsns = c()

for (t in 1:timepoints) {

	# update of probabilities and state transitions

 	if (egg==1) {p1=1} else {p1=0}

 	if (vc==0 && runif(1)<p1) {vcn = 1}      	               # vc switches from 0 to 1
 	if (vc==1 && runif(1)<p2) {vcn = 0}           	         # vc switches from 1 to 0
 
 	if (count>=clusterN) {p3=1} else {p3=0}

 	if (uv1==0 && runif(1)<p3) {uv1n = 1} 	        	   # uv1 change from 0 to 1
 	if (uv1==1 && runif(1)<p4) {uv1n = 0}                    # uv1 change from 1 to 0

 	if ((uv1+vc)==0) {p5=1} else {p5=0}

 	if (runif(1)<p5) {hsnn=1} else {hsnn=0}
 	if (hsnn==1) {p6=1} else {p6=0}                       

 	if (egg==0 && runif(1)<p6) {eggn=1}			         # egg-laying
 
 	if (egg==1 && runif(1)<p7) {eggn=0}

 	count <- p0*count + eggn

	# update state of the system
 	egg = eggn; vc = vcn; uv1 = uv1n; hsn = hsnn  

	# output state variables for new time point
 	eggs <- cbind(eggs,egg)
 	vcs <- cbind(vcs,vc)
 	uv1s <- cbind(uv1s,uv1) 
 	counts <- cbind(counts, count)
 	hsns<-cbind(hsns,hsn)
};

# plot results
plot(c(1:timepoints)*tau, eggs, type="l", col="black", lwd=2, xlab="time/[sec]", ylab="egg laying events")
lines(c(1:timepoints)*tau, vcs/2, col="red", lwd=2)
lines(c(1:timepoints)*tau, uv1s/3, col="green", lwd=2)

# calculate egg-laying statistics
eggnumber[k]=sum(eggs)/timepoints/tau;           # the egg-laying rate (per second)
b=(which(uv1s==0))
c=(which(uv1s[b+1]==0))
b=b[-c]
CN[k]=sum(eggs[b[1]:b[length(b)]])/(length(b)-1) # the average number of egg-laying events per cluster
}

