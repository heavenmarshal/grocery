# Code for STAT 464/864 Final Exam - Solutions
#Francois Marshall, Queen's University
##################################################################################

#This code is to be used as an example of the standard which is required in any assignment which requires programming.
#Commenting is an important part of programming.  It is used to make rigorous what it is that you are trying to do, not only for
#recollection upon future readings of the code for the user, but also for anyone else who may need to read it.


library(Hmisc)                          #used in q2
library(multitaper)
library(astsa)
library(locfit)
library(zoo)
library(xts)
library(outliers)
library(e1071)
library(FBN)
library(Rmpfr)
library(IDPmisc)
library(date)
library(moments)
library(phonTools)
library(wq)
library(gtools)
library(formattable)


#Set working directory, so that the plot will not be saved under the home directory; rather, it should be saved in the current, or working directory.
setwd("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Final_Exam/Solutions")



###############################################################################
###############################################################################
#User-defined functions.

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











###############################################################################
###############################################################################




#Where the output goes.

sink("Final_Exam_Solutions_Output.txt", append=FALSE, split=FALSE)




##################################################################################
##################################################################################
# Question 1
##################################################################################
##################################################################################

cat("Question 1\n\n\n")


dataset<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Yearly_Data.txt", header=TRUE)

temperature.series<-dataset$Temperature
temperature.series<-ts(temperature.series)

N.samples=length(temperature.series)

initial.date=1764
sampling.period=1/12 #Years.
dates<-1:N.samples
dates<-dates*sampling.period+initial.date



#Remove the multitaper estimate of the mean value.

NW=3
K=2*NW-1
M=2^(log2(N.samples)+4)

cat("NW = ",NW,", K = ",K,", M = ",M,"\n\n")

#Compute the tapers.
slepian.sequences<-dpss.matrix(N.samples,NW,K)
#Compute the Slepian functions.
slepian.functions<-Slepian.functions(slepian.sequences,M)
#Compute the eigencoefficients.
eigencoeffs<-eigencoefficients(temperature.series,slepian.sequences,M)

#De-trend the estimate of the mean term for the time series.
mean.estimate=multitaper.mean_estimate(eigencoeffs,slepian.functions)
temperature.series<-temperature.series-mean.estimate






##################################################################################
##################################################################################
# Part b)
##################################################################################
##################################################################################


#Read in the significant frequencies.

frequency.dataset<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Final_Exam/Significant_Frequencies.txt", header=TRUE)

significant_frequencies<-frequency.dataset$Frequency
significant.F_statistics<-frequency.dataset$F_statistic
significant.periods<-frequency.dataset$Period

n.significant=length(significant_frequencies)

#To do: Sort both the frequencies and F-statistics in order of increasing F-statistic.  Include the output in a table.

#FM: This is the solution.

ordered.indices<-rev(order(significant.F_statistics))

orderedFrequencies<-formattable(significant_frequencies[ordered.indices],digits=3,format="f")
orderedFStatistics<-formattable(significant.F_statistics[ordered.indices],digits=3,format="f")
orderedPeriods<-formattable(significant.periods[ordered.indices],digits=3,format="f")

frequencies.table<-data.frame(orderedFrequencies,orderedFStatistics,orderedPeriods)
write.csv(frequencies.table,file = "Significant_Frequencies.csv",row.names = FALSE)










##################################################################################
##################################################################################
# Part c)
##################################################################################
##################################################################################


#The following code is used to read in the necessary numbers for computation of the complex-demodulate reconstruction
#of the 1-year periodic component.

#As stated in the complex-demodulate reconstruction formula in Thomson 1995 (the equation in the sentence immediately
#after Equation (3b)), the demodulates must be multiplied by complex exponentials.  These complex exponentials are
#included in the "year.complex_exponential" vectors.




##################################################################################

#Read in the complex demodulates.

real.matrix<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Complex_Demodulates_Real.txt", header=TRUE)
year.re_demodulate<-real.matrix$X7
#To do: Read in the sixth and thirteenth columns.
#FM: This is the solution.
year.re_demodulate1<-real.matrix$X13
year.re_demodulate2<-real.matrix$X6

imaginary.matrix<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Complex_Demodulates_Imaginary.txt", header=TRUE)
year.im_demodulate<-imaginary.matrix$X7
#To do: Read in the sixth and thirteenth columns.
#FM: This is the solution.
year.im_demodulate1<-imaginary.matrix$X13
year.im_demodulate2<-imaginary.matrix$X6

year.complex_demodulate<-year.re_demodulate+1i*year.im_demodulate
#To do: Repeat for the other two demodulates.
#FM: This is the solution.
year.complex_demodulate1<-year.re_demodulate1+1i*year.im_demodulate1
year.complex_demodulate2<-year.re_demodulate2+1i*year.im_demodulate2



#Read in the complex exponentials.

real.matrix<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Complex_Exponentials_Real.txt", header=TRUE)
year.re_exponential<-real.matrix$X7
#To do: Read in the sixth and thirteenth columns.
#FM: This is the solution.
year.re_exponential1<-real.matrix$X13
year.re_exponential2<-real.matrix$X6

imaginary.matrix<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Complex_Exponentials_Imaginary.txt", header=TRUE)
year.im_exponential<-imaginary.matrix$X7
#To do: Read in the sixth and thirteenth columns.
#FM: This is the solution.
year.im_exponential1<-imaginary.matrix$X13
year.im_exponential2<-imaginary.matrix$X6

year.complex_exponential<-year.re_exponential+1i*year.im_exponential
#To do: Repeat for the other two complex exponentials.
#FM: This is the solution.
year.complex_exponential1<-year.re_exponential1+1i*year.im_exponential1
year.complex_exponential2<-year.re_exponential2+1i*year.im_exponential2

##################################################################################




#Fourier-series amplitudes, used to weigh the demodulates.

mean.values_object<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Mean_Values.txt", header=TRUE)

mean.values_re<-mean.values_object$Real
mean.values_im<-mean.values_object$Imaginary
complex.mean_values<-complex(real=mean.values_re,imaginary=mean.values_im)
mean_values<-Mod(complex.mean_values)





##################################################################################

#Compute a complex-demodulate reconstruction.

N.demodulate=length(year.complex_demodulate)
block.length=12*5 #Five years, as recommended in Thomson, 1995.
demodulate.indices<-1:(N.demodulate-block.length) #Needed two blocks per one-block demodulate, so ignore the last block.

#From Thomson, 1995.
multiplied.sequences<-year.complex_exponential[demodulate.indices]*year.complex_demodulate[demodulate.indices]
reconstructed.cosine=2*Re(multiplied.sequences)
#Need to scale the cosine functions so that, together, they add up to well approximate the time series.
reconstruction<-mean_values[7]*reconstructed.cosine #Not in Thomson, 1995, but shows the importance of the wave in the superposition.

#To do: Repeat for the other two complex exponentials.
#FM: This is the solution.
multiplied.sequences1<-year.complex_exponential1[demodulate.indices]*year.complex_demodulate1[demodulate.indices]
reconstructed.cosine1<-2*Re(multiplied.sequences1)
reconstruction1<-mean_values[13]*reconstructed.cosine1

multiplied.sequences2<-year.complex_exponential2[demodulate.indices]*year.complex_demodulate2[demodulate.indices]
reconstructed.cosine2<-2*Re(multiplied.sequences2)
reconstruction2<-mean_values[6]*reconstructed.cosine2


##################################################################################




#Plot the complex-demodulate reconstructions, overlaid on the original time series.

truncated.temperature_series<-temperature.series[demodulate.indices]

min.temperature=min(reconstruction)
max.temperature=max(reconstruction)

x_plotting<-dates[demodulate.indices]
y_plotting<-truncated.temperature_series

lower_x=1820
upper_x=1840
lower_y=min(truncated.temperature_series[x_plotting>=lower_x & x_plotting<=upper_x])
upper_y=22

pdf("Demodulate_Reconstruction.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Date, in years",ylab="Temperature, in degrees Celsius",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting,y_plotting,lwd=2)
lines(x_plotting,reconstruction,col=2,lty=1,lwd=2)
#To do: Repeat for the other two complex exponentials.
#FM: This is the solution.
lines(x_plotting,reconstruction1,col=3,lty=1,lwd=2)
lines(x_plotting,reconstruction2,col=4,lty=1,lwd=2)
leg.text<-c("Temperature series","0.27-year cycle","1-year cycle","1.18-year cycle")
legend(1832,22,leg.text,lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c(1,3,2,4))
dev.off()












##################################################################################
##################################################################################
# Question 2
##################################################################################
##################################################################################

cat("Question 2\n\n\n")

#Read in the exchange-rate series.

dataset<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Final_Exam/Exchange_Rates_Full.txt", header=TRUE)

dates<-dataset$Time
exchange.rates<-dataset$Rate

N.exchange=length(exchange.rates)




##################################################################################
##################################################################################
# Part a)
##################################################################################
##################################################################################

#To do: Determine and plot the time series of the first and second differences of the
#original exchange-rate time series.

#FM: This is the solution.

x_plotting<-dates
y_plotting<-exchange.rates

lower_x=min(x_plotting)
upper_x=max(x_plotting)
lower_y=min(y_plotting)
upper_y=max(y_plotting)

pdf("Exchange_Rates.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Days since April 4th, 2016",ylab="US dollar to Canadian dollar daily exchange rate",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting,y_plotting,lwd=2)
dev.off()


#Compute the first differences.
first.differences<-exchange.rates[2:N.exchange]-exchange.rates[1:(N.exchange-1)]
N.first=length(first.differences)

#Compute the second differences.
second.differences<-first.differences[2:N.first]-first.differences[1:(N.first-1)]
N.second=length(second.differences)


#Plot the second differences.

x_plotting<-dates[1:N.second]
y_plotting<-second.differences

lower_x=min(x_plotting)
upper_x=max(x_plotting)
lower_y=min(y_plotting)
upper_y=max(y_plotting)

pdf("Second_Differences.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Days since April 4th, 2016",ylab="Second difference of exchange rate, in USD/CAD",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting,y_plotting,lwd=2)
dev.off()





#Read in the table of first and second differences which is given.
dataset<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Final_Exam/Differenced_Series.txt", header=TRUE)

threshold=-1

second.differences<-dataset$Second
second.differences<-second.differences[second.differences>threshold]

N.second=length(second.differences)





##################################################################################
##################################################################################
# Part b)
##################################################################################
##################################################################################

NW=5
K=2*NW-1
M=2^(log2(N.second)+2)

#Compute the tapers.
slepian.sequences<-dpss.matrix(N.second,NW,K)
#Compute the Slepian functions.
slepian.functions<-Slepian.functions(slepian.sequences,M)
#Compute the eigencoefficients.
eigencoeffs<-eigencoefficients(second.differences,slepian.sequences,M)


#Compute the frequencies.
frequency.obj<-frequencies.principle_domain(M)
frequencies<-frequency.obj$temp.frequencies1
frequencies<-non_negative.principle_domain(frequencies)

#Compute the F-statistic.
indices<-1:M
mu.values<-sapply(indices,harmonic.mean,slepian.par=slepian.functions,eigencoefficent.par=eigencoeffs)
mu.powers<-sapply(indices,harmonic.signal_power,slepian.par=slepian.functions,mu.par=mu.values)
mu.residuals<-sapply(indices,harmonic.residual,slepian.par=slepian.functions,eigencoefficent.par=eigencoeffs,mu.values)
F.statistics<-harmonic.F_test(mu.powers,mu.residuals)
F.statistics<-non_negative.principle_domain(F.statistics)

##################################################################################



#Compute the significant frequencies.
M.2=length(F.statistics)
threshold.quantile=qf(0.95,2,2*K-2)
peak.frequencies<-c()
peak.heights<-c()
peak.indices<-c()
counter=1
for(i in 2:(M.2-1)){
  if(F.statistics[i]>threshold.quantile & F.statistics[i]>F.statistics[i-1] & F.statistics[i]>F.statistics[i+1]){
    peak.frequencies[counter]=frequencies[i]
    peak.heights[counter]=F.statistics[i]
    peak.indices[counter]=i
    counter=counter+1
  }
}

n.significant=length(peak.frequencies)

##################################################################################


#Compute the complex-demodulate reconstructions of the narrow-band functions centred on
#the significant frequencies.

mu.values<-non_negative.principle_domain(mu.values)
N.2=N.second/2
x.deterministic<-rep(0,N.2)
time.indices<-1:N.second-1
time_indices.half<-1:N.2-1

for (i in 1:n.significant) {
  temp.exponentials<-exp(-1i*2*pi*peak.frequencies[i]*time.indices)
  temp.signal<-temp.exponentials*as.complex(second.differences)
  temp.filter<-as.complex(slepian.sequences[,1])
  temp.transform<-fft(temp.signal)*fft(temp.filter)
  temp.demodulate<-fft(temp.transform,inverse=TRUE)
  temp.demodulate<-temp.demodulate[1:N.2]
  temp.cosine<-2*Re(exp(1i*2*pi*peak.frequencies[i]*time_indices.half)*temp.demodulate)
  temp.index=peak.indices[i]
  temp.reconstruction<-Mod(mu.values[temp.index])*temp.cosine/(2*NW)
  x.deterministic<-x.deterministic+temp.reconstruction
}




##################################################################################

#Only use the first half of the series, since that is the domain which the complex demodulate covers.
dates<-dates[1:N.2]
second.differences<-second.differences[1:N.2]
N.second=N.2

#Remove the periodic trend from the series.
stationary.series<-second.differences-x.deterministic


##################################################################################

#Plot the sequence of partial autocorrelation coefficients.

pdf("Partial_Autocorrelations.pdf",width=8,height=6)
pacf(stationary.series)
dev.off()



##################################################################################

#Compute a multitaper estimate of the spectrum.

M=2^(log2(N.second)+4)
M.2=M/2

#Compute the frequencies.
frequency.obj<-frequencies.principle_domain(M)
frequencies<-frequency.obj$temp.frequencies1
frequencies<-non_negative.principle_domain(frequencies)

#Compute the tapers.
slepian.sequences<-dpss.matrix(N.second,NW,K)
#Compute the Slepian functions.
slepian.functions<-Slepian.functions(slepian.sequences,M)
#Compute the eigencoefficients.
eigencoeffs<-eigencoefficients(stationary.series,slepian.sequences,M)
#Compute the eigenspectra.
eigenspectra<-Mod(eigencoeffs)^2

#Compute the multitaper estimate of the spectrum.
multitaper.spectrum<-multitaper.average_spectrum(eigenspectra)
#Take only the non-negative half of the estimate.
multitaper_spectrum.non_negative<-non_negative.principle_domain(multitaper.spectrum)
multitaper.spectrum<-c(multitaper_spectrum.non_negative,rev(multitaper_spectrum.non_negative))








##################################################################################
##################################################################################
# Part c)
##################################################################################
##################################################################################

#Compute a multitaper estimate of the autocorrelation sequence.

autocovariance.sequence<-Re(fft(multitaper.spectrum,inverse=TRUE))/M
autocorrelation.sequence<-autocovariance.sequence/autocovariance.sequence[1]
autocorrelation.sequence<-autocorrelation.sequence[1:M.2]


#To do: Compute and plot the AR-2 fit for the residual time series.

#FM: This is the solution.

a_11=autocorrelation.sequence[2]
a_12=autocorrelation.sequence[2]*(1-autocorrelation.sequence[3])/(1-autocorrelation.sequence[2]^2)
a_22=(autocorrelation.sequence[3]-autocorrelation.sequence[2]^2)/(1-autocorrelation.sequence[2]^2)

cat("Estimates for the first two partial-correlation coefficients are ",a_12," and ",a_22,".\n\n")

sample.var=mean(multitaper.spectrum)
var_2=(1-a_11^2)*(1-a_22^2)*sample.var

cat("The estimate for the 2-step prediction variance is ",var_2,".\n\n")

var_2.iid=(1-a_11^2)*(1-a_22^2)*var(stationary.series)

cat("Using the standard sample variance estimate, ",var(stationary.series),", for an independent and independently-distributed sample, the estimate for the 2-step prediction variance is ",var_2.iid,".\n\n")



ar2.fit<-sapply(1:N.second,ar2.point_estimate,ts.par=stationary.series,coeff.par=c(a_12,a_22))
ar2.fit<-ar2.fit[3:N.second]


#Not required in the solution, but useful to include dispersion characteristics in the exam solution set.

acov_matrix<-toeplitz.matrix(autocovariance.sequence[1:2])
inverse.acov_matrix<-solve(acov_matrix)
scaled.iacov_matrix<-var_2*inverse.acov_matrix

phi.parameters<-c(a_11,a_12)
phi.LB<-phi.parameters-qnorm(0.95)*sqrt(diag(scaled.iacov_matrix)/N.second)
phi.UB<-phi.parameters+qnorm(0.95)*sqrt(diag(scaled.iacov_matrix)/N.second)

ar2.LB<-sapply(1:N.second,ar2.point_estimate,ts.par=stationary.series,coeff.par=phi.LB)
ar2.LB<-ar2.LB[3:N.second]
ar2.UB<-sapply(1:N.second,ar2.point_estimate,ts.par=stationary.series,coeff.par=phi.UB)
ar2.UB<-ar2.UB[3:N.second]


x_plotting<-dates
y_plotting<-stationary.series

lower_x=min(x_plotting)
upper_x=max(x_plotting)
lower_y=min(y_plotting)
upper_y=2*max(y_plotting)

pdf("AR2_Fit_Dispersion.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Days since April 4th, 2016",ylab="Exchange rate, in USD/CAD (gross trend components removed)",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting,y_plotting,type="o",pch=".",cex=4,lwd=3)
lines(x_plotting[3:N.second],ar2.LB,col=4,lwd=2)
lines(x_plotting[3:N.second],ar2.UB,col=7,lwd=2)
lines(x_plotting[3:N.second],ar2.fit,col=2,lwd=2)
leg.text<-c("Exchange rates","AR-2 fit","AR-2 lower-bound parameters","AR-2 upper-bound parameters")
legend(0,0.07,leg.text,fill=c(1,2,4,7))
dev.off()


ar2.LB<-ar2.fit-qnorm(0.95)*sqrt(var_2)
ar2.UB<-ar2.fit+qnorm(0.95)*sqrt(var_2)


x_plotting<-dates
y_plotting<-stationary.series

lower_x=min(x_plotting)
upper_x=max(x_plotting)
lower_y=min(y_plotting)
upper_y=2*max(y_plotting)

pdf("AR2_Fit.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Days since April 4th, 2016",ylab="Exchange rate, in USD/CAS (gross trend components removed)",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting[3:N.second],ar2.LB,pch=".",col="grey")
lines(x_plotting[3:N.second],ar2.UB,pch=".",col="grey")
polygon(c(x_plotting[3:N.second],rev(x_plotting[3:N.second])),c(ar2.UB,rev(ar2.LB)),col="grey75",border=NA)
lines(x_plotting,y_plotting,type="o",pch=".",cex=4,lwd=3)
lines(x_plotting[3:N.second],ar2.fit,col=2,lwd=2)
leg.text<-c("Exchange rates","AR-2 fit","95% confidence interval")
legend(0,0.07,leg.text,fill=c(1,2,"grey75"))
dev.off()



##################################################################################
##################################################################################
# Part e)
##################################################################################
##################################################################################

#Read in the fitted series.
dataset<-read.table("/Users/francoismarshall/Documents/Frank's folder/Queen's/TA/2016-17/Winter/STAT 464_864/webpage/Assignments/Final/Final_Exam/AR2_Fit.txt", header=TRUE)
ar2.fit<-dataset$Fit
num.predictions=length(ar2.fit)


#To do: Subtract the AR-2 fitted series, "ar2.fit", from the residual series, "stationary.series[3:N.second]"
#Compute the multitaper estimate of the spectrum.

#FM: This is the solution.

#Compute the spectrum for the nearly-white series.

stationary.series<-stationary.series[3:N.second]-ar2.fit

N.residual=length(stationary.series)

M=2^(log2(N.residual)+4)
M.2=M/2

#Compute the frequencies.
frequency.obj<-frequencies.principle_domain(M)
frequencies<-frequency.obj$temp.frequencies1
frequencies<-non_negative.principle_domain(frequencies)

#Compute the tapers.
slepian.sequences<-dpss.matrix(N.residual,NW,K)
#Compute the Slepian functions.
slepian.functions<-Slepian.functions(slepian.sequences,M)
#Compute the eigencoefficients.
eigencoeffs<-eigencoefficients(stationary.series,slepian.sequences,M)
#Compute the eigenspectra.
eigenspectra<-Mod(eigencoeffs)^2

#Compute the multitaper estimate of the spectrum.
multitaper.spectrum<-multitaper.average_spectrum(eigenspectra)
#Take only the non-negative half of the estimate.
multitaper_spectrum.non_negative<-non_negative.principle_domain(multitaper.spectrum)

x_plotting<-frequencies
y_plotting<-log(multitaper_spectrum.non_negative)

lower_x=min(x_plotting)
upper_x=max(x_plotting)
lower_y=-14#min(y_plotting)
upper_y=5#max(y_plotting)

pdf("Multitaper_Spectrum_White.pdf",width=8,height=6)
par(mgp=c(2,1,0))
plot(0,0,xlab="Frequency, in cycles per day",ylab="Log power spectral density, (spectrum in squared USD/CAD times days)",main="",xlim=c(lower_x,upper_x),ylim=c(lower_y,upper_y),pch=".",lab=c(10,10,7))
grid()
minor.tick(10,10)
lines(x_plotting,y_plotting,lwd=2)
dev.off()









sink()


