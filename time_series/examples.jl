# Create some standard time series and plot them.
using Plots: plot, savefig

nsamples = 1000
normal_samples = randn(nsamples)
pushfirst!(normal_samples, 0.0)

# White Noise
white_noise = cumsum(normal_samples)
fig = plot(1:length(white_noise), white_noise, label="White Noise")
savefig(fig, "white_noise.png")

# Moving Average

