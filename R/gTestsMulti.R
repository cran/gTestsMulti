## main function
gtestsmulti = function(E, data_list, perm=0){
  K = length(data_list)
  n = sapply(data_list, nrow)
  N = sum(n)
  ind = 1:N

  R = getR(E, n, ind)

  Ebynode = vector("list", N)
  for (i in 1:N) Ebynode[[i]] = rep(0,0)
  for (i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]], E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]], E[i,1])
  }
  nE = nrow(E)
  nodedeg = rep(0,N)
  for (i in 1:N) nodedeg[i] = length(Ebynode[[i]])
  nEi = sum(nodedeg*(nodedeg-1))  # pair of nodes sharing a node * 2

  muo = rep(0, K)
  sdo = rep(0, K)

  p1 = n*(n-1)/N/(N-1)
  p2 = p1*(n-2)/(N-2)
  p3 = p2*(n-3)/(N-3)

  quan = nE*(nE-1)-nEi
  muo = nE*p1
  sdo = sqrt( nEi*p2 + quan*p3 + muo - muo^2 )

  mu = diag(K)
  cov11 = diag(K) # covariance function of R_{ii}
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      mu[i,j] = nE*2*n[i]*n[j]/N/(N-1)
      cov11[i,j] = ( quan*p1[i]*n[j]*(n[j]-1)/(N-2)/(N-3) - nE^2*p1[i]*p1[j] )
    }
  }
  cov11[lower.tri(cov11)] = t(cov11)[lower.tri(cov11)]
  diag(cov11) = sdo^2
  diag(mu) = muo

  R[lower.tri(R)] = t(R)[lower.tri(R)] # edge matrix
  mu[lower.tri(mu)] = t(mu)[lower.tri(mu)] # mean matrix

  # covariance of R_{ii} and R_{jk}
  num = K*(K-1)/2
  cov12 = matrix(0, K, num)
  temp = t(combn(K,2))
  for (i in 1:K) {
    temp1 = matrix(i, num, 2)
    temp2 = sapply(1:num, function(x) setdiff(temp[x,],temp1[x,]))
    for (j in 1:num) {
      temp3 = temp2[[j]]
      if (length(temp3)==1) {
        cov12[i,j] = nEi*n[i]*(n[i]-1)*n[temp3]/N/(N-1)/(N-2) + 2*quan*n[i]*(n[i]-1)*(n[i]-2)*n[temp3]/N/(N-1)/(N-2)/(N-3) - mu[i,i]*mu[i,temp3]
      } else {
        cov12[i,j] = quan*2*n[i]*(n[i]-1)*prod(n[temp3])/N/(N-1)/(N-2)/(N-3) - mu[i,i]*mu[temp3[1],temp3[2]]
      }
    }
  }

  # variance of R_{ij}
  sdo1 = rep(0, num)
  for (i in 1:num) {
    temp1 = n[temp[i,]]
    temp2 = prod(temp1)
    temp3 = sum(temp1)
    temp4 = mu[lower.tri(mu)]^2
    sdo1[i] = nE*2*temp2/N/(N-1) + nEi*temp2*(temp3-2)/N/(N-1)/(N-2) + quan*4*temp2*(temp2-temp3+1)/N/(N-1)/(N-2)/(N-3) - temp4[i]
  }

  # covariance of R_{ij} and R_{kl}
  temp1 = apply(temp, 1, prod)
  temp2 = apply(temp, 1, sum)
  temp3 = temp1 - temp2 + 1
  cov22 = matrix(0, num, num)
  for (i in 1:num) {
    u = temp[i,]
    for (j in i:num) {
      v = temp[j,]
      temp1 = intersect(u,v)
      temp2 = union(u,v)
      if (length(temp2)==3) {
        cov22[i,j] = nEi*prod(n[temp2])/N/(N-1)/(N-2) + quan*4*prod(n[temp2])*(n[temp1]-1)/N/(N-1)/(N-2)/(N-3) - mu[u[1],u[2]]*mu[v[1],v[2]]
      } else {
        cov22[i,j] = quan*4*prod(n[temp2])/N/(N-1)/(N-2)/(N-3) - mu[u[1],u[2]]*mu[v[1],v[2]]
      }
    }
  }
  cov22[lower.tri(cov22)] = t(cov22)[lower.tri(cov22)]
  diag(cov22) = sdo1

  Sig = rbind(cbind(cov11, cov12), cbind(matrix(0, K*(K-1)/2, K), cov22))
  Sig[lower.tri(Sig)] = t(Sig)[lower.tri(Sig)] # covariance matrix of whole R_{ii} and R_{ij}
  KK = nrow(Sig)

  inv = solve(cov11)
  Sw = ( diag(R)-diag(mu) ) %*% inv %*% ( diag(R)-diag(mu) )
  S_W = Sw[1] # SW
  
  if ( class(try(solve(cov22),silent=T))[1] != "matrix" ) {
      inv1 = ginv(cov22)
  } else {
      inv1 = solve(cov22)
  }
  vec1 = R[lower.tri(R)] - mu[lower.tri(mu)]
  Sb = vec1 %*% inv1 %*% vec1
  S_B = Sb[1] # SB

  S = S_W + S_B # S

  cov33 = Sig[1:(KK-1),1:(KK-1)]
  if ( class(try(solve(cov33),silent=T))[1] != "matrix" ) {
      inv2 = ginv(cov33)
  } else {
      inv2 = solve(cov33)
  }
  vec2 = c(diag(R),R[lower.tri(R)][-K*(K-1)/2]) - c(diag(mu),mu[lower.tri(mu)][-K*(K-1)/2])
  Sa = vec2 %*% inv2 %*% vec2
  S_A = Sa[1] # SA

  S_W_pval.approx = pchisq(S_W, df=K, lower.tail=F)
  S_B_pval.approx = pchisq(S_B, df=rankMatrix(cov22), lower.tail=F)
  S_pval.approx = 2*min( S_W_pval.approx, S_B_pval.approx )
  S_A_pval.approx = pchisq(S_A, df=rankMatrix(cov33), lower.tail=F)


  result = list()
  result$teststat$S = S # test statistic S
  result$teststat$S_A = S_A # test statistic S_A
  result$pval$S_appr = min(1,S_pval.approx) # approximated p-value of S
  result$pval$S_A_appr = min(1,S_A_pval.approx) # approximated p-value of S_A

  if (perm > 0) {
    S_W_perm = S_B_perm = S_perm = S_A_perm = rep(0, perm)
    for (i in 1:perm) {
      sam = sample(N)
      R_perm = getR(E, n, sam)

      R_perm[lower.tri(R_perm)] = t(R_perm)[lower.tri(R_perm)]

      Sw_perm = ( diag(R_perm)-diag(mu) ) %*% inv %*% ( diag(R_perm)-diag(mu) )
      S_W_perm[i] = Sw_perm[1]

      vec1_perm = R_perm[lower.tri(R_perm)] - mu[lower.tri(mu)]
      Sb_perm = vec1_perm %*% inv1 %*% vec1_perm
      S_B_perm[i] = Sb_perm[1]

      S_perm[i] = S_W_perm[i] + S_B_perm[i]

      vec2_perm = c(diag(R_perm),R_perm[lower.tri(R_perm)][-K*(K-1)/2]) - c(diag(mu),mu[lower.tri(mu)][-K*(K-1)/2])
      Sa_perm = vec2_perm %*% inv2 %*% vec2_perm
      S_A_perm[i] = Sa_perm[1]
    }
    S_pval.perm = length(which(S_perm>=S))/perm
    S_A_pval.perm = length(which(S_A_perm>=S_A))/perm

    result$pval$S_perm = min(1,S_pval.perm) # permutation p-value of S
    result$pval$S_A_perm = min(1,S_A_pval.perm) # permutation p-value of S_A
  }

  return( result )
}


## supporting function for constructing R
getR = function(E, n, ind){
  K = length(n)
  R = matrix(0, K, K)
  E_ind = rep(1:K, n)
  for (i in 1:nrow(E)) {
    e1 = E_ind[match(E[i,1], ind)]
    e2 = E_ind[match(E[i,2], ind)]
    R[e1, e2] = R[e1, e2] + 1
  }
  R[upper.tri(R)] = t(R)[upper.tri(R)] + R[upper.tri(R)]
  R[lower.tri(R)] = 0
  return(R)
}
