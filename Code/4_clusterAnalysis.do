version 17.0

set seed 1234

egen sin24_std = std(sin24_hat)
egen cos24_std = std(cos24_hat)
egen sin12_std = std(sin12_hat)
egen cos12_std = std(cos12_hat)
egen sin8_std  = std(sin8_hat)
egen cos8_std  = std(cos8_hat)
egen mesor_std  = std(mesor_hat)

cluster kmeans sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std, k(7) gen(kGroup) measure(L2)

drop *_std