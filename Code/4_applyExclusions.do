version 17.0

// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05
drop if totalPAEE_p >0.05

drop if mesor_hat == .
drop if amplitude24_hat == . | amplitude12_hat == . | amplitude8_hat == .
drop if acrophase24_hat == . | acrophase12_hat == . | acrophase8_hat  == .

drop if acrophase24_se == .
drop if acrophase12_se == .
drop if acrophase8_se == .

drop if P_Pwear_consolidated<72
drop if fatMass == .
drop if fatFreeMass == .


