Benchmark CTDirect for different AD backends (keyword AD_backend = ...)
- :default : ForwardDiff
- :optimized (default for CTDirect) : Forward / ReverseDiff
- :enzyme : Enzyme
- :zygote : Zygote

Using default benchmark

Grid size       default     optimized       enzyme      zygote
250
500
1000
2500
5000