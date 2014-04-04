#plot expected abundance for the different data generating models

Covariate=c(0:100)/100
Exp.N.RD=exp(3+10*Covariate-10*Covariate^2)  #restricted dispersal - relationship

pdf('Covariate_abundance_relation.pdf')
plot(Covariate,Exp.N.RD,type="l",lwd=2,xlab="Covariate",ylab="Abundance relationship")
dev.off()