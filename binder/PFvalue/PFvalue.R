library("markovchain")

#Matriz de interés A para determinar su eigenvalor dominante (Perron-Frobenius).
#Cambiarla según las necesidades.

A = matrix(data = c(10, 8, 2,
                    0.3, 4, 1,
                    2, 15, 5), byrow = T, nrow = 3)

eigen(A)
###########################################################

#Cadena de Markov asociada a la matriz de interés A

tam = dim(A)[1]
sumA = rowSums(A)
maxA = max(sumA)

transAux = cbind(A/maxA,1-sumA/maxA)
transAux = rbind(transAux,c(rep(0,tam),1))

States = c(paste("s",1:tam,sep=""),"ext")

mcAux <- new("markovchain", states = States,
             transitionMatrix = transAux, name = "Aux")
###########################################################

#Simulación Monte Carlo

N = 1e3
mcSize = 500 #tamaño de muestra -cadenas de Markov-
z = States[1] #seleccionar el estado retorno
tau = rep(0,N)
tauS2 = rep(0,N)
tauS3 = rep(0,N)
varS2 = rep(0,N)
varS3 = rep(0,N)

for(i in 1:N)
{
 simul <- rmarkovchain(n = mcSize, object = mcAux, t0 = z)
 idx1 = which(simul==z)
 idx2 = which(simul=="ext")
 if(length(idx1)!=0 && idx1[1]<idx2[1]) tau[i]=idx1[1]; 

 simS2 <- rmarkovchain(n = mcSize, object = mcAux, t0 = States[2])
 simS3 <- rmarkovchain(n = mcSize, object = mcAux, t0 = States[3])
 
 idx1 = which(simS2==z)
 idx2 = which(simS2=="ext")
 if(length(idx1)!=0 && idx1[1]<idx2[1]) tauS2[i]=idx1[1]; 

 idx1 = which(simS3==z)
 idx2 = which(simS2=="ext")
 if(length(idx1)!=0 && idx1[1]<idx2[1]) tauS3[i]=idx1[1]; 
} #fin for i in 1:N

#Bajo la aproximación sum(exp(n_k*x)) = exp(N_max*x). n={1,2,3,..,} y N_max = max{n_k} 
log(N)/(range(tau)[2])

#Solución de la ecuación sum(exp(tau_k*x)) = N
y = function(x){sum(exp(tau*x)) - N}
b = unlist(lapply(0:10,function(x) y(x))) #Evalúa la función y en [0,10]
b = which(b>0)[1] #Detecta el 1er valor x en [0,10] tal que y(x)>0
theta = uniroot(y, interval = c(0, b))
 
uS2 = sum(exp(tauS2*theta$root))/N
uS3 = sum(exp(tauS3*theta$root))/N

#varianza muestral de uS2 y uS3
varS2 = (sum(exp(2*tauS2*theta$root))/N - (uS2)^2)*(1/(N-1)) 
varS3 = (sum(exp(2*tauS3*theta$root))/N - (uS3)^2)*(1/(N-1))
 
v = c(1,uS2,uS3)

#valor estimado del eigenvalor dominante de la matriz A y un múltiplo de su vector asociado
PFvector = as.vector(A%*%v)
PFvalue = PFvector[1]

#error relativo de uS2 y uS3
if(uS2>0) REs2 <- sqrt(varS2)/uS2 else REs3<- NA 
if(uS3>0) REs3 <- sqrt(varS3)/uS3 else REs3<- NA 
###########################################################

#Definición de varias funciones

f <- function(u)
{
  transIS = mcAux[-(tam+1),-(tam+1)]
  transIS = matrix(mapply(function(x,y) x*y,t(transIS),u), byrow=T, nrow=tam)
  transIS = transIS/rowSums(transIS)
  transIS = cbind(transIS,rep(0,tam))
  transIS = rbind(transIS,c(rep(0,tam),1))
 
  mc <- new("markovchain", states = States,
            transitionMatrix = transIS, name = "IS")
  mc
}

LH <- function(x,z,MC,AUX,IS)
{
  aux = c(x,MC)
  aux1 = which(aux==z)
  if(x==z) aux = aux[1:(aux1[2]-1)] else aux = aux[1:(aux1[1]-1)]
  MC1 = MC[1:which(MC==z)[1]]
  PTrans = unlist(mapply(function(x,y) AUX[x,y],aux,MC1))
  QTrans = unlist(mapply(function(x,y) IS[x,y],aux,MC1))
  names(PTrans) = NULL
  names(QTrans) = NULL
 
  if(length(MC1)>0)
  {
    tau = length(MC1)
    LHRatio = prod(PTrans/QTrans)
    output <- list("tau" = tau, "LHRatio" = LHRatio)
  } else output <-list("tau" = 0, "LHRatio" = 0)
  output
}

#Cross-Entropy             
LH_CE <- function(x,z,MC,AUX,IS)
{
  aux=which(MC==z)
  SimAux=MC[1:aux[1]] #truncar la cadena hasta el primer regreso a z
  aux=which(SimAux==x)
  if(length(aux) == 0) n1=0 else 
  {
    SimAux1=SimAux[aux] #visitas al estado x
    n1=length(SimAux1)
  }
  
  SimAux1 = SimAux[aux+1]
  rem = which(is.na(SimAux1)==T)
  if(length(rem) > 0) SimAux1 = SimAux1[-rem]
  if(length(SimAux1)==0) n2=rep(list(0),tam+1) else
  {
    SimAux2=lapply(States,function(x) which(SimAux1==x)) #visitas a cada estado
    n2=lapply(SimAux2,function(x) length(x))
  }
  
  aux = c(z,SimAux)
  aux1 = which(aux==z)
  aux = aux[1:(aux1[2]-1)]
  if(length(aux1)==0) MC1=NULL else 
  {
    MC1 = SimAux
    PTrans = unlist(mapply(function(x,y) AUX[x,y],aux,MC1))
    QTrans = unlist(mapply(function(x,y) IS[x,y],aux,MC1))
    names(PTrans) = NULL
    names(QTrans) = NULL
  }
  
  if(length(MC1) > 0)
  {
    LHRatio = prod(PTrans/QTrans)
    out1 <- n1*LHRatio
    out2 <- unlist(lapply(n2,function(x) x*LHRatio))
    output <- list("Est1" = out1, "Est2" = out2)
  } else output <-list("Est1" = 0,"Est2" = rep(0,tam+1))
  output
}

g <- function(st,List,txt)
{  
  if(txt == "tau") res = unlist(lapply(List,function(x) x[[st]]$tau))
  if(txt == "LHRatio") res = unlist(lapply(List,function(x) x[[st]]$LHRatio))
  res
}
#######################################################################################

#Simulación Monte Carlo vía IS adaptado de acuerdo al algoritmo de Desai, Glynn (2001).

N = 1e2
M = N/10
mcSize = 500 #tamaño de muestra -cadenas de Markov-
z = States[1] #seleccionar el estado retorno
v = rep(1,tam) #vector inicial
PFvalue = rep(0,N)
lst3 = list("tau"=0,"LHRatio"=0)
lst2 = rep(list(lst3),tam)
names(lst2) = States[-(tam+1)]
lst = rep(list(lst2),M)

for(i in 1:N)
{
  for(k in 1:M)
  {
      mcIS = f(v) #Cadena de Markov con matriz de transisión Q_v
      Sim <- lapply(States[-(tam+1)],function(x) rmarkovchain(n = mcSize, object = mcIS, t0 = x))
      Out <- mapply(function(x1,x2) LH(x1,z,x2,mcAux,mcIS),States[-(tam+1)],Sim)

      #Aplicar la función por columnas a la matriz Out
      lst[[k]] = apply(Out,2,function(x) x)
  }#fin for k in 1:M
 
  PosRetSt = which(States==z)
  tau = lapply(States[-(tam+1)],function(x) g(x,lst,"tau"))
  LHRatio = lapply(States[-(tam+1)],function(x) g(x,lst,"LHRatio"))
 
  #Solución (x) de la ecuación sum(LHRatio_k*exp(tau_k*x)) = M iniciando en el estado retorno z
  y = function(x){sum(LHRatio[[PosRetSt]]*exp(tau[[PosRetSt]]*x)) - M}
  b = unlist(lapply(0:10,function(x) y(x))) #Evalúa la función y en [0,10]
  b = which(b>0)[1] #Detecta el 1er valor x en [0,10] tal que y(x)>0
  theta = uniroot(y, interval = c(0, b))
 
  #Actualización del vector v
  v = mapply(function(x,y) sum(x*exp(y*theta$root))/M,LHRatio,tau)
  #v[PosRetSt]=1

  #Estimaciones parciales del eigenvalor dominante de A 
  #PFvalue[i] = as.numeric(A[PosRetSt,]%*%v)
  PFvalue[i] = maxA*exp(-theta$root)

}#fin for i in 1:N

#Estimación del eigenvalor dominante de A 
PFvalue[N]
############################################################################

#Simulación Monte Carlo vía IS con Cross-Entropy
             
M = 1e4
N = M/10
mcSize = 500 #tamaño de muestra -cadenas de Markov-
z = States[1] #seleccionar el estado retorno

#v = rep(1,tam) #vector inicial
#mcIS = f(v) #Cadena de Markov con matriz de transisión Q_v

#Matriz de transición inicial para IS
transIS <- matrix(data = rep(1/tam,tam*tam), byrow = T, nrow = 3)
transIS <- cbind(transIS,rep(0,tam))
transIS <- rbind(transIS,c(rep(0,tam),1))
mcIS <- new("markovchain", states = States,
                      transitionMatrix = transIS, name = "IS")

alpha = 0.6
Est1 = rep(list(0),tam+1)
Est2 = rep(list(c(0,0,0,0)),tam+1)

lst3 = list("tau"=0,"LHRatio"=0)
lst2 = rep(list(lst3),tam)
names(lst2) = States[-(tam+1)]
lst = rep(list(lst2),N)

#Cross-Entropy. Iteraciones para aproximar la matriz de transición óptima
for(i in 1:M)
{ 
  z0 = sample(States[-(tam+1)],1)
  Sim <- rmarkovchain(n = mcSize, object = mcIS, t0 = z0)
  Out <- lapply(States, function(x) LH_CE(x,z0,Sim,mcAux,mcIS))
  
  a = lapply(Out,function(x) x$Est1)
  b = lapply(Out,function(x) x$Est2)
  Est1 = as.list(mapply(function(x,y) x+y,Est1,a))
  Est2 = mapply(function(x,y) x+y,Est2,b)
  Est2 = lapply(seq_len(ncol(Est2)), function(i) Est2[,i])
  mcISAux = lapply(seq_len(nrow(mcIS[,])), function(i) mcIS[i,])
  transIS = alpha*t(mapply(function(x,y,z) if(y > 0 && sum(x/y) == 1) x/y else z,Est2,Est1,mcISAux)) + (1-alpha)*mcIS[,] 
  mcIS <- new("markovchain", states = States,
              transitionMatrix = transIS, name = "IS")
}

#Iteraciones para estimar el eigenvalor dominante de la matriz A                           
for(i in 1:N)
{
  Sim <- lapply(States[-(tam+1)],function(x) rmarkovchain(n = mcSize, object = mcIS, t0 = x))
  Out <- mapply(function(x1,x2) LH(x1,z,x2,mcAux,mcIS),States[-(tam+1)],Sim)
    
  #Aplicar la función por columnas a la matriz Out
  lst[[i]] = apply(Out,2,function(x) x)
}#fin for i in 1:N

  PosRetSt = which(States==z)
  tau = lapply(States[-(tam+1)],function(x) g(x,lst,"tau"))
  LHRatio = lapply(States[-(tam+1)],function(x) g(x,lst,"LHRatio"))
  
  #Solución (x) de la ecuación sum(LHRatio_k*exp(tau_k*x)) = N iniciando en el estado retorno z
  y = function(x){sum(LHRatio[[PosRetSt]]*exp(tau[[PosRetSt]]*x)) - N}
  b = unlist(lapply(0:10,function(x) y(x))) #Evalúa la función y en [0,10]
  b = which(b>0)[1] #Detecta el 1er valor x en [0,10] tal que y(x)>0
  theta = uniroot(y, interval = c(0, b))
  
  #Estimaciones parciales del eigenvalor dominante de A 
  PFvalue = maxA*exp(-theta$root)
  
#Estimación del eigenvalor dominante de A 
PFvalue
############################################################################
