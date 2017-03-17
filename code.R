## Simulate the sample paths
mu=.01
sigma=.2
lambda=2
k=0
s=.1
m=5
n=250
dt=1/n
S=matrix(log(20000),m,n+1)
for (i in 1:n) {
	jump=ifelse(runif(m)<lambda*dt,1,0)
	jumpsize=jump*rnorm(m,k,s)
	S[,i+1]=S[,i]+rnorm(m,mu*dt,sigma*sqrt(dt))+jumpsize
}
plot(exp(S[1,]),type='l',xlab='time', ylab='equity price', ylim=c(min(exp(S)),max(exp(S))))
for (i in 2:m) {
	points(exp(S[i,]), type='l')
}

## Gibbs sampler
Y=diff(log(S[2,]))
# Set the initial values
mu=0
sigma=.1
lambda=1
k=0.1
s=.5
jump=c(rep(1,n))
jumpsize=Y/2
# Assign the prior
mean.mu=mean(Y)/dt/n
var.mu=1
alpha.sigma=1
beta.sigma=.01
alpha.lambda=1
beta.lambda=.01
mean.k=.5
var.k=.5
alpha.s=1
beta.s=.01

M=10000
est.mu=0
est.sigma=0
est.lambda=0
est.k=0
est.s=0
for (i in 1:M){
	# Calculate posterior
	v.mu=1/(1/(sigma^2/n/dt)+1/var.mu)
	m.mu=((sum(Y-jumpsize)/n/dt)/(sigma^2/n/dt)+mean.mu/var.mu)/(1/(sigma^2/n/dt)+1/var.mu)
	mu=rnorm(1,m.mu,sqrt(v.mu))
	
	a.sigma=n/2+alpha.sigma
	b.sigma=beta.sigma+sum((Y-mu*dt-jumpsize)^2)/2
	sigma=sqrt(rinvgamma(1,a.sigma,b.sigma))/sqrt(dt)
	
	J=jumpsize[jump==1]
	j=sum(jump)
	a.lambda=j+alpha.lambda
	b.lambda=n-j+beta.lambda
	lambda=rbeta(1,a.lambda,b.lambda)/dt
	
	if (j>1){
		v.k=1/(1/(s^2/j)+1/var.k)
		m.k=(mean(J)/(s^2/j)+mean.k/var.k)/(1/(s^2/j)+1/var.k)
		k=rnorm(1,m.k,sqrt(v.k))
		
		a.s=j/2+alpha.s
		b.s=beta.s+sum((J-k)^2)/2
		s=sqrt(rinvgamma(1,a.s,b.s))
	}
	
	pjump=1/(1+(1-lambda*dt)/(lambda*dt)*sqrt((sigma^2*dt+s^2)/(sigma^2*dt))*exp(-(Y-mu*dt)^2/(sigma^2*dt)/2+(Y-mu*dt-k)^2/(sigma^2*dt+s^2)/2))
	jump=ifelse(runif(n)<pjump,1,0)
	var.jump=1/(1/s^2+1/sigma^2/dt)
	mean.jump=((Y-mu*dt)/sigma^2/dt+k/s^2)*var.jump
	jumpsize=jump*(rnorm(n)*sqrt(var.jump)+mean.jump)
	est.mu=est.mu+mu/M
	est.sigma=est.sigma+sigma/M
	est.lambda=est.lambda+lambda/M
	est.k=est.k+k/M
	est.s=est.s+s/M
}
c(est.mu,est.sigma,est.lambda,est.k,est.s)

## Hang Seng Index
dat=read.csv("HSI.csv")
HSI=rev(dat[,5])
plot.ts(HSI)
Y=diff(log(HSI))
n=length(Y)
dt=1/n
# Then use the same code as in the simulation case
