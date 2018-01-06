frequencies.principle_domain<-function(M.par,sampling_period.par=1,N.par=1){
  #Returns the frequencies for a spectrum with M samples between the Nyquist bounds.
  indices<-1:M.par
  temp.frequencies<-indices-M.par/2 #Centre the frequencies.
  temp.frequencies<-temp.frequencies/M.par #Scale the frequencies to {-1/2+1/M,...,1/2}
  temp.frequencies1<-temp.frequencies/sampling_period.par #Scale the frequencies to lie between the Nyquist bounds
  #(e.g., set to units of cycles per year).
  temp.frequencies2<-temp.frequencies*N.par #Multiplying by N to convert to Rayleighs.
  frequencies.object<-data.frame(temp.frequencies1,temp.frequencies2)
  return(frequencies.object)
}

non_negative.principle_domain<-function(vector.par){
  #Return the values of the vector which are associated with non-negative frequencies.
  temp.M=length(vector.par)
  temp.M_2=temp.M/2
  return(vector.par[temp.M_2:temp.M])
}

return.to.fft_ordering<-function(vector.par){
  #Return the vector to the order of the FFT output.
  temp.M=length(vector.par)
  temp.M_2=temp.M/2
  temp.vector<-rep(0,temp.M)
  temp.vector[1:(temp.M_2+1)]<-vector.par[temp.M_2:temp.M]
  temp.vector[(temp.M_2+2):temp.M]<-vector.par[1:(temp.M_2-1)]
  return(temp.vector)
}

pad.matrix<-function(matrix.par,M.par){
  print(M.par)
  temp.N=nrow(matrix.par)
  temp.K=ncol(matrix.par)
  temp.matrix<-matrix.par
  num.padded=M.par-temp.N
  temp.matrix<-rbind(temp.matrix,matrix(0,num.padded,temp.K))
  return(temp.matrix)
}

fft.rearrange_matrix<-function(vector.par){
  #Rearrange the FFT output to have the same ordering as in Blackman and Tukey: {-1/2+1/M,...,1/2}.
  temp.M=nrow(vector.par)
  temp.K=ncol(vector.par)
  temp.M_2=temp.M/2
  temp.matrix<-matrix(0,nrow=temp.M,ncol=temp.K)
  temp.matrix[1:(temp.M_2-1),]<-vector.par[(temp.M_2+2):temp.M,]
  temp.matrix[temp.M_2:temp.M,]<-vector.par[1:(temp.M_2+1),]
  return(temp.matrix)
}

dpss.matrix<-function(N.par,NW.par,K.par){
  cat("dpss.matrix: N.par = ",N.par,", NW.par = ",NW.par,", K.par = ",K.par,"\n\n")
  dpss.object<-dpss(N.par,K.par,NW.par)
  dpss.tapers<-dpss.object$v
  return(dpss.tapers)
}

Slepian.functions<-function(dpss.par,M.par){
  temp.N=nrow(dpss.par)
  temp.K=ncol(dpss.par)
  padded.matrix<-pad.matrix(dpss.par,M.par)
  cat("M = ",nrow(padded.matrix),"\n\n")
  slepian.fft<-mvfft(padded.matrix)
  slepian.fft<-fft.rearrange_matrix(slepian.fft)
  return(slepian.fft)
}

Slepian.windows<-function(slepian.par){
  slepian.windows<-Mod(slepian.par)^2
  #slepian.windows[!slepian.windows]<-max(slepian.windows)*1e-12
  return(slepian.windows)
}

eigencoefficients<-function(ts.par,dpss.par,M.par){
  temp.N=nrow(dpss.par)
  temp.K=ncol(dpss.par)
  tapered.series<-t(t(dpss.par) %*% diag(ts.par))
  padded.matrix<-pad.matrix(tapered.series,M.par)
  eigencoefficient.fft<-mvfft(padded.matrix)
  eigencoefficient.fft<-fft.rearrange_matrix(eigencoefficient.fft)
  return(eigencoefficient.fft)
}

multitaper.average_spectrum<-function(eigenspectra.par){
  return(rowMeans(eigenspectra.par))
}

delete.one_matrix<-function(k.par,matrix.par){
  temp.num_rows=nrow(matrix.par)
  temp.num_columns=ncol(matrix.par)
  temp.matrix<-matrix(0,nrow=temp.num_rows,ncol=temp.num_columns-1)
  if(k.par==1){
    temp.matrix<-matrix.par[,2:temp.num_columns]
  }
  else if(k.par==temp.num_columns){
    temp.matrix<-matrix.par[,1:(temp.num_columns-1)]
  }
  else{
    temp.matrix<-cbind(matrix.par[,1:(k.par-1)],matrix.par[,(k.par+1):temp.num_columns])
  }
  return(temp.matrix)
}

delete.one_spectrum<-function(k.par,eigenspectra.par){
  temp.matrix<-delete.one_matrix(k.par,eigenspectra.par)
  temp.spectrum<-multitaper.average_spectrum(temp.matrix)
  return(temp.spectrum)
}

jackknife.variance<-function(matrix.par){
  temp.num_rows=nrow(matrix.par)
  temp.num_columns=ncol(matrix.par)
  delete.one_spectra<-matrix(0,nrow=temp.num_rows,ncol=temp.num_columns)
  for(k in 1:temp.num_columns){
    delete.one_spectra[,k]<-delete.one_spectrum(k,matrix.par)
  }
  temp.log_spectra<-log(delete.one_spectra)
  temp.log_averages<-rowMeans(temp.log_spectra)
  temp.difference_matrix<-sweep(temp.log_spectra,MARGIN=1,temp.log_averages,FUN="-")
  temp.log_variances<-(temp.num_columns-1)*rowMeans(temp.difference_matrix^2)
  temp.dt<-data.frame(temp.log_averages,temp.log_variances)
  return(temp.dt)
}

jackknife.bounds<-function(matrix.par,alpha.par){
  temp.num_columns=ncol(matrix.par)
  temp.jackknife_object<-jackknife.variance(matrix.par)
  temp.quantile=qt(1-alpha.par/2,temp.num_columns-1)
  temp.log_avgs<-temp.jackknife_object$temp.log_averages
  temp.log_vars<-temp.jackknife_object$temp.log_variances
  temp.log_sigmas<-sqrt(temp.log_vars)
  temp.log_LB<-temp.log_avgs-temp.quantile*temp.log_sigmas
  temp.log_UB<-temp.log_avgs+temp.quantile*temp.log_sigmas
  temp.LB<-exp(temp.log_LB)
  temp.UB<-exp(temp.log_UB)
  temp.dt<-data.frame(temp.log_avgs,temp.LB,temp.UB,temp.log_LB,temp.log_UB)
  return(temp.dt)
}

multitaper.acvs<-function(spectrum.par){
  #Computes the multitaper estimate of the autocovariance function.
  #Returns both the autocovariance sequence and the associated lags, where lags are in units of the sampling period.
  temp.M=length(spectrum.par)
  temp.M_2=temp.M/2
  temp.spectrum<-return.to.fft_ordering(spectrum.par)
  temp.acvs<-fft(temp.spectrum,inverse=TRUE)/temp.M
  temp.acvs<-Re(temp.acvs)
  temp.acvs<-temp.acvs[1:(temp.M_2+1)]
  temp.lags<-1:M-1
  temp.lags<-temp.lags[1:(temp.M_2+1)]
  temp.acvs_object<-data.frame(temp.acvs,temp.lags)
  return(temp.acvs_object)
}

multitaper.acrs<-function(spectrum.par){
  temp.multitaper_acvs_object<-multitaper.acvs(spectrum.par)
  temp.acvs<-temp.multitaper_acvs_object$temp.acvs
  temp.lags<-temp.multitaper_acvs_object$temp.lags
  temp.variance<-temp.acvs[1]
  temp.acrs<-temp.acvs/temp.variance
  temp.acrs_object<-data.frame(temp.acvs,temp.acrs,temp.lags)
  return(temp.acrs_object)
}


toeplitz.matrix<-function(acvs.par){
  n.lags=length(acvs.par)
  temp.matrix<-matrix(0,nrow=n.lags,ncol=n.lags)
  for(j in 1:(n.lags-1)){
    temp.vector<-c(acvs.par[j:1],acvs.par[2:(n.lags-j+1)])
    #cat(j,", ",length(temp.vector),", ",n.lags,"\n")
    temp.matrix[,j]<-temp.vector
  }
  temp.matrix[,n.lags]<-rev(acvs.par)
  return(temp.matrix)
}


multitaper.mean_estimate<-function(eigencoefficients.par,Slepian_function.par){
  temp.M=nrow(eigencoefficients.par)
  zero.index=temp.M/2
  temp.K=ncol(eigencoefficients.par)
  even.indices<-seq(from=1,to=temp.K,by=2)
  y_k.0<-eigencoefficients.par[zero.index,]
  y_k.0<-y_k.0[even.indices]
  U_k.0<-Re(Slepian_function.par[zero.index,])
  U_k.0<-U_k.0[even.indices]
  temp.numerator_sum=crossprod(y_k.0,U_k.0)
  temp.denominator.sum=crossprod(U_k.0)
  temp.mean_estimate=Re(temp.numerator_sum/temp.denominator.sum)
  return(temp.mean_estimate)
}

harmonic.mean<-function(frequency.par,slepian.par,eigencoefficent.par){
  temp.M=nrow(slepian.par)
  temp.zero_index=temp.M/2
  temp.numerator<-sum(slepian.par[temp.zero_index,]*eigencoefficent.par[frequency.par,])
  temp.denominator<-sum(Mod(slepian.par[temp.zero_index,])^2)
  temp.mu<-temp.numerator/temp.denominator
  return(temp.mu)
}

harmonic.signal_power<-function(frequency.par,slepian.par,mu.par){
  temp.M=nrow(slepian.par)
  temp.zero_index=temp.M/2
  temp.mu=mu.par[frequency.par]
  temp.mod_mu<-Mod(slepian.par[temp.zero_index,])
  temp.signal_power<-Mod(temp.mu)^2*crossprod(temp.mod_mu)
  return(temp.signal_power)
}

harmonic.residual<-function(frequency.par,slepian.par,eigencoefficent.par,mu.par){
  temp.M=nrow(slepian.par)
  temp.zero_index=temp.M/2
  temp.differences<-eigencoefficent.par[frequency.par,]-mu.par[frequency.par]*slepian.par[temp.zero_index,]
  temp.squared_differences<-Mod(temp.differences)^2
  temp.rss<-sum(temp.squared_differences)
  return(temp.rss)
}

harmonic.F_test<-function(signal.par,residual.par){
  temp.F=(K-1)*signal.par/residual.par
  return(temp.F)
}

running.sum<-function(end.par,sequence.par){
  temp.sum=0
  if(end.par==1){
    temp.sum=sequence.par[1]
  }
  else{
    temp.sum=sum(sequence.par[1:end.par])
  }
  return(temp.sum)
}

serial.convolution<-function(sequence1_par,sequence2_par){
  temp.N=length(sequence1_par)
  temp.indices1<-1:(2*temp.N-1)
  temp.indices2<-sapply(temp.indices1,convolution.index,start_par=-temp.N)
  temp.convolution<-sapply(temp.indices2,multiply.sequences,sequence1_par=sequence1_par,sequence2_par=sequence2_par)
  return(temp.convolution)
}

complex.demodulate<-function(ts.par,freq.par,filter.par){
  temp.N=length(ts.par)
  temp.indices<-1:temp.N-1
  temp.harmonic_arguments<-2*pi*freq.par*temp.indices
  temp.cosines<-cos(temp.harmonic_arguments)
  temp.sines<-sin(temp.harmonic_arguments)*(-1)
  re_ts.par<-ts.par*temp.cosines
  im_ts.par<-ts.par*temp.sines
  complex.series<-complex(real=re_ts.par,imaginary=im_ts.par)
  temp.fft<-fft(complex.series)
  temp.filter_transform<-fft(filter.par)
  transformed.fft<-temp.fft*temp.filter_transform
  temp.demodulate<-fft(transformed.fft,inverse=TRUE)/temp.N
  temp.N_2=temp.N/2
  temp.demodulate<-temp.demodulate[1:temp.N_2]
  temp.times<-temp.indices[1:temp.N_2]
  temp.demodulate_object<-data.frame(temp.times,temp.demodulate)
  return(temp.demodulate_object)
}

harmonic.peaks<-function(frequencies.par,F.par,K.par,threshold_percentile.par=0.99){
  threshold_quantile.par=10*log10(qf(threshold_percentile.par,2,2*K.par-2))
  temp.N=length(F.par)
  temp.decibel_F<-10*log10(F.par)
  #Find the maxima using the cumulative sum.
  index.sequence<-1:temp.N
  log.sums<-sapply(index.sequence,running.sum,temp.decibel_F)
  temp.first_diffs<-log.sums[2:temp.N]-log.sums[1:(temp.N-1)]
  temp.N_first<-length(temp.first_diffs)
  temp.second_diffs<-temp.first_diffs[2:temp.N_first]-temp.first_diffs[1:(temp.N_first-1)]
  temp.threshold_diff=0
  where.positive<-which(temp.second_diffs>temp.threshold_diff)
  indices.prior_inflection<-list()
  counter=1
  num.optimal=length(where.positive)
  for(i in 1:num.optimal){
    temp.index=where.positive[i]
    if(temp.index<temp.N & temp.second_diffs[temp.index+1]<=temp.threshold_diff){
      indices.prior_inflection[[counter]]=temp.index+2
      counter=counter+1
    }
  }
  optimal.indices<-as.numeric(indices.prior_inflection)
  num.prior_inflection=length(optimal.indices)
  optimal.frequencies<-frequencies.par[optimal.indices]
  maxima<-temp.decibel_F[optimal.indices]
  peak.minima_indices<-matrix(0,nrow=num.prior_inflection,ncol=2)
  for(i in 1:num.prior_inflection){
    left.sample<-c()
    right.sample<-c()
    temp.centre_index=0
    temp.left_index=0
    temp.right_index=0
    if(i==1){
      temp.left_index=optimal.indices[i]
      temp.right_index=optimal.indices[i+1]
      left.sample<-temp.decibel_F[1:(temp.left_index-1)]
      right.sample<-temp.decibel_F[(temp.left_index+1):(temp.right_index-1)]
    }
    else if(i==num.prior_inflection){
      temp.left_index=optimal.indices[i-1]
      temp.right_index=optimal.indices[i]
      left.sample<-temp.decibel_F[(temp.left_index+1):(temp.right_index-1)]
      right.sample<-temp.decibel_F[(temp.right_index+1):temp.N]
    }
    else{
      temp.left_index=optimal.indices[i-1]
      temp.centre_index=optimal.indices[i]
      temp.right_index=optimal.indices[i+1]
      left.sample<-temp.decibel_F[(temp.left_index+1):(temp.centre_index-1)]
      right.sample<-temp.decibel_F[(temp.centre_index+1):(temp.right_index-1)]
    }
    temp.left_minimum=min(left.sample)
    temp.right_minimum=min(right.sample)
    if(i==1){
      temp.left_minimum_index=max(which(left.sample==temp.left_minimum))
    }
    else{
      temp.left_minimum_index=temp.left_index+max(which(left.sample==temp.left_minimum))
    }
    temp.right_minimum_index=temp.centre_index+min(which(right.sample==temp.right_minimum))
    peak.minima_indices[i,1]=temp.left_minimum_index
    peak.minima_indices[i,2]=temp.right_minimum_index
  }
  minimum.frequencies<-frequencies.par[peak.minima_indices[,1]]
  minima<-temp.decibel_F[peak.minima_indices[,1]]
  minimum.frequencies2<-frequencies.par[peak.minima_indices[,2]]
  minima2<-temp.decibel_F[peak.minima_indices[,2]]
  optimal.freq<-rep(0,num.prior_inflection)
  optimal.F_statistics<-rep(0,num.prior_inflection)
  for(i in 1:num.prior_inflection){
    temp.centre_index=optimal.indices[i]
    temp.LB=peak.minima_indices[i,1]
    temp.UB=peak.minima_indices[i,2]
    peak.frequencies<-frequencies.par[temp.LB:temp.UB]
    peak.F_statistics<-temp.decibel_F[temp.LB:temp.UB]
    temp.minimum=0
    if(temp.decibel_F[temp.LB]<temp.decibel_F[temp.UB]){
      temp.minimum=temp.decibel_F[temp.UB]
    }
    else{
      temp.minimum=temp.decibel_F[temp.LB]
    }
    peak.frequencies<-peak.frequencies[peak.F_statistics>=temp.minimum]
    peak.F_statistics<-peak.F_statistics[peak.F_statistics>=temp.minimum]
    #Compute the centre of mass.
    temp.com=crossprod(peak.frequencies,peak.F_statistics)/sum(peak.F_statistics)
    temp.com_F_statistic=max(peak.F_statistics)
    frequency.at_max=peak.frequencies[peak.F_statistics==max(peak.F_statistics)]
    right.middle_frequency=median(peak.frequencies[peak.frequencies>=min(peak.frequencies) & peak.frequencies<=frequency.at_max])
    left.middle_frequency=median(peak.frequencies[peak.frequencies>=frequency.at_max & peak.frequencies<=max(peak.frequencies)])
    if(temp.com>=right.middle_frequency || temp.com<=left.middle_frequency){
      optimal.freq[i]=crossprod(peak.frequencies,peak.F_statistics)/sum(peak.F_statistics)
      optimal.F_statistics[i]=max(peak.F_statistics)
    }
    else{
      optimal.freq[i]=frequency.at_max
      optimal.F_statistics[i]=peak.F_statistics[peak.F_statistics==max(peak.F_statistics)]
    }
  }
  temp.significant_frequencies<-optimal.freq[optimal.F_statistics>threshold_quantile.par]
  temp.significant_F_statistics<-optimal.F_statistics[optimal.F_statistics>threshold_quantile.par]
  output.object<-data.frame(temp.significant_frequencies,temp.significant_F_statistics)
  return(output.object)
}

ar2.point_estimate<-function(ts.par,coeff.par,index.par){
  temp.N=length(ts.par)
  temp.estimate=0
  if(index.par>2){
    temp.estimate=coeff.par[1]*ts.par[index.par-1]+coeff.par[2]*ts.par[index.par-2]
  }
  else{
    temp.estimate=0
  }
  return(temp.estimate)
}

ar.denominator<-function(freq.par,coeff.par){
  temp.order=length(coeff.par)
  exp.vector<-complex(real=cos(2*pi*(0:temp.order)*freq.par),imaginary=-sin(2*pi*(0:temp.order)*freq.par))
  temp.out<-Mod(sum(c(1,-coeff.par)*exp.vector))^2
  return(temp.out)
}

ar.spectrum<-function(var.par,N.par,coeff.par){
  #Large M for good resolution.
  temp.M=2^(log2(N.par)+2)
  temp.frequencies<-1:(temp.M/2+1)/temp.M-1
  temp.spectrum<-var.par/sapply(temp.frequencies,ar.denominator,coeff.par=coeff.par)
  return(temp.spectrum)
}

harmonic.reconstruction<-function(index.par,frequencies.par,amplitudes.par){
  temp.amplitudes<-Mod(amplitudes.par)
  temp.phases<-atan2(Im(amplitudes.par),Re(amplitudes.par))
  temp.cosines<-temp.amplitudes*cos(2*pi*frequencies.par*index.par+temp.phases)
  temp.sum<-2*sum(temp.cosines)
  return(temp.sum)
}
