la1 = 2 
la2 = 3 
gam = 3
psi = 0.5

r1 <- c(la1^-1+la1^(gam-2)*psi*gam*(2*gam-1),0,gam*psi*la1^(gam-1),gam*la1^(gam-1))
r2 <- c(0, la2^-1+la2^(gam-2)*psi*gam*(2*gam-1),gam*psi*la2^(gam-1),gam*la2^(gam-1))
r3 <- c(gam*psi*la1^(gam-1),gam*psi*la2^(gam-1),psi*(la1^gam+la2^gam),la1^gam+la2^gam)
r4 <- c(gam*la1^(gam-1),gam*la2^(gam-1),la1^gam+la2^gam,(la1^gam+la2^gam)/psi)

mat <- matrix(c(r1,r2,r3,r4),nrow=4)

det(mat)



la1 = 2 
la2 = 3 
la3 = 1
gam = 3
psi = 0.5

r1 <- c(la1^-1+la1^(gam-2)*psi*gam*(2*gam-1),0,0,gam*psi*la1^(gam-1),gam*la1^(gam-1))
r2 <- c(0, la2^-1+la2^(gam-2)*psi*gam*(2*gam-1),0,gam*psi*la2^(gam-1),gam*la2^(gam-1))
r3 <- c(0, 0,la3^-1+la3^(gam-2)*psi*gam*(2*gam-1),gam*psi*la3^(gam-1),gam*la3^(gam-1))
r4 <- c(gam*psi*la1^(gam-1),gam*psi*la2^(gam-1),gam*psi*la3^(gam-1),psi*(la1^gam+la2^gam+la3^gam),la1^gam+la2^gam+la3^gam)
r5 <- c(gam*la1^(gam-1),gam*la2^(gam-1),gam*la3^(gam-1),la1^gam+la2^gam+la3^gam,(la1^gam+la2^gam+la3^gam)/psi)


mat <- matrix(c(r1,r2,r3,r4,r5),nrow=5)

det(mat)
