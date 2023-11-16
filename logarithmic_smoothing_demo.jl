### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ c86650f9-4952-4064-80b6-e6d4c6253126
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
end

# ╔═╡ dfb41db2-0b57-4de9-a44e-cd1f70842588
using Plots, CSV, ImageFiltering, NonuniformResampling1D, NonlinearSequences

# ╔═╡ e2ce0c30-0ce1-11ed-084b-4523774d8541
md"""
# Logarithmic Smoothing in Julia

This demo shows how to smooth data on a logarithmic scale in Julia. As an example, I'll use a loudspeaker frequency response measurement.

The smoothing will use two unregistered packages of mine: [`NonuniformResampling1D`](https://github.com/Firionus/NonuniformResampling1D.jl) and [`NonlinearSequences`](https://github.com/Firionus/NonlinearSequences.jl).

The HTML version of this notebook can be viewed at <https://firionus.github.io/logarithmic_smoothing_julia_demo>. 

To run this Pluto notebook yourself, see the instructions at <https://github.com/Firionus/logarithmic_smoothing_julia_demo>.

Let's get started!
"""

# ╔═╡ ddfb9f83-a845-42d9-801a-d866528c3603
md"""
## Import the Packages

For this notebook, the supplied `Project.toml` is used. 

To perform logarithmic smoothing in your own environment, install the two unregistered packages. For this, switch your REPL into Pkg mode by typing `]` and run:
```
add https://github.com/Firionus/NonuniformResampling1D.jl
add https://github.com/Firionus/NonlinearSequences.jl
```
"""

# ╔═╡ cc004309-0ed9-465d-9f43-037ab1949a4e
md"""
## Read the Example Data

The example file is an in-room frequency response measurement of a 2.5" midrange/broadband loudspeaker (Tymphany PLS-65F25AL02-08 in approximately 2l closed box). It was exported as CSV from the acoustic measurement software [ARTA](https://artalabs.hr/). 

Let's use CSV.jl to parse it, convert each column to an array and save the frequencies as `f` and the magnitudes in dB as `magnitude`:
"""

# ╔═╡ d8a3e8f1-83b4-425c-ace8-c48452c97c73
f, magnitude, _ = map( # ignore phase in third column
	# exclude first value (f=0) since it's no good on a log scale 
	col -> (col.column|>Array)[2:end], 
	CSV.File("example_frequency_response.csv").columns
)

# ╔═╡ 2d0aba3c-b2e2-4ea1-8912-9e0d75532080
md"""
## Plot the Raw Data

The magnitude response data for loudspeakers is typically plotted in dB with a logarithmic frequency scale from the lowest audible frequency of 20 Hz to the highest audible frequency of 20 kHz. The logarithmic scale is important because it approximates the way humans perceive tones. Fundamental frequencies of 220 Hz, 440 Hz and 880 Hz are commonly heard as the same tone, just with an octave between each of them (a factor of two). 

> Comment: There are better scales for approximating the human perception of tones, like mel, EBR or Bark. But none of them are as simple as the logarithmic scale. 

Looking at this plot, the data is quite clear around 150 Hz. The output drops off below 120 Hz, which means the speaker does not produce much bass, which is expected from a midrange driver. 

However, the interpretation of the high frequencies is significantly hindered by the fine grained noise at high frequencies above 1 kHz. The German loudspeaker legend Anselm Goertz likes to call it "grass". Because of the logarithmic x-axis, this problem only appears at the high frequencies, where the data points on the linear frequency scale are compressed together. 

Another issue is the large size of the measurement data, which contains $(length(f)) points. This resolution could be much reduced at high frequencies. 
"""

# ╔═╡ 1963ee15-c177-4e64-aa22-60fb783a1ec3
plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Magnitude Response - Raw Data", label="raw", xlims=(20, 20e3), ylims=(40, 120))

# ╔═╡ df40186e-d36a-4779-be83-f709e102c43d
md"""
### Failing with Linear Smoothing

Let's use a naive moving mean smoothing algorithm from [`ImageFiltering.jl`](https://juliaimages.org/ImageFiltering.jl/stable/) and see what happens: 
"""

# ╔═╡ 35d2be48-81b5-4605-882a-684deb8675d3
begin
	plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Linear Smoothing Doesn't Work", label="raw", xlims=(20, 20e3), ylims=(40, 120), alpha=.3)
	plot!(f, imfilter(magnitude, KernelFactors.box(199)), label="lin smooth $(round(199*f[1])|>Int) Hz width")
	plot!(f, imfilter(magnitude, KernelFactors.box(9_999)), label="lin smooth $(round(9_999*f[1])|>Int) Hz width")
end

# ╔═╡ a4dfa7f8-b727-4936-bb24-63521168535a
md"""
This shows the issue with linear smoothing: Either the noise at high frequencies is gone, but the low frequencies have no detail and have issues with missing data outside the boundaries (that's why the magnitude goes down below 1 kHz). Or the low frequencies have the right amount of detail, but the highs are too noisy. 

Further, the filtering for the larger filter takes over a second, which is unacceptably slow. This is because there are so many data points to consider. Therefore, the smoothing algorithm should be combined with downsampling to reduce memory requirements and speed up computation. 

That is where logarithmic smoothing and resampling comes in. 

## Logarithmic Smoothing

First, we'll define our desired resolution in octaves:
"""

# ╔═╡ a44c0978-a733-4a90-b871-0c844586dcec
resolution = 1/12 # octaves

# ╔═╡ 80fab818-246f-4c69-bd43-faa23d4b3a16
md"""
Next, define the desired frequency points after the smoothing and resampling. We use [`NonlinearSequences.jl`](https://github.com/Firionus/NonlinearSequences.jl) to create a logarithmically spaced frequency vector:
"""

# ╔═╡ b255464e-366e-4329-b6c0-4d619d2be4ff
flog = octspace(20, 21e3, resolution)

# ╔═╡ 9f6810e4-eadb-410b-850f-6bb9cda78800
md"""
The resampling requires the input values to be sampled linearly and the input frequencies to be defined as a range. This makes calculations efficient and restricts the resampling problem to a manageable subset.  So let's approximate the linearly spaced frequencies with a range:
"""

# ╔═╡ 11f8da69-4a42-4f5c-a77c-11d94181c0a8
frange = range(f[1], f[end], length=length(f))

# ╔═╡ 5b1b6194-53f6-42b7-ac77-c691a373df3b
md"""
Now we use `nuresample` from [`NonuniformResampling1D.jl`](https://github.com/Firionus/NonuniformResampling1D.jl) to resample the high resolution input data to our desired `flog` points. The bandwidth for the smoothing is automatically chosen based on the distance between the points in `flog` and is therefore logarithmic. 
"""

# ╔═╡ dbfc2f92-d52a-4b04-b554-8d4c491d7549
smoothed = nuresample(frange, magnitude, flog)

# ╔═╡ 6def680a-5cff-474f-a556-eb311f118dba
begin
	plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Logarithmic Smoothing", label="raw", xlims=(20, 20e3), ylims=(40, 120), alpha=.3)
	plot!(flog, smoothed, label="log smooth 1/12 oct")
end

# ╔═╡ 7956290b-3a19-40c2-8194-ff444d2b0e6a
md"""
Now this is a much nicer result! Both the low and the high frequencies are smoothed but retain the same amount of detail when plotted on a log axis. 

And because we performed downsampling at the same time, `smoothed` is only $(length(smoothed)) elements long and is calculated in less than a millisecond. Over 1000 times faster than linear smoothing!

Of course, I hear you say, the smoothed graph now looks quite edgy because the amount of points is so low. Luckily, there's an easy way to change that: oversampling.
"""

# ╔═╡ 8a332d89-b63a-43a2-af02-907cfbd11505
md"""
## Making It Look Nice with Oversampling

Oversampling increases the amount of points without changing the width of the smoothing. Let's define the oversampling factor:
"""

# ╔═╡ db47c349-4dc6-431a-8081-286506a674e3
oversampling = 8 # 8 times as many points as previously

# ╔═╡ ebf1d7f5-ca3b-47fb-8fab-e49104a4fae5
flog_oversampled = octspace(20, 21e3, resolution/oversampling)

# ╔═╡ 6d4e2046-4d64-4267-81a6-29df79763542
md"""
Increase thewidth of the `smoothing_function` by the oversampling factor:
"""

# ╔═╡ 7b56eec9-2717-49df-b39b-017b4d7a330c
oversampled = nuresample(frange, magnitude, flog_oversampled, rect_window(.5*oversampling))

# ╔═╡ f67b1de1-1247-4ecf-be84-0bc1816bad34
begin
	plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Logarithmic Smoothing with Oversampling", label="raw", xlims=(20, 20e3), ylims=(40, 120), alpha=.3)
	plot!(flog_oversampled, oversampled, label="log smooth 1/12 oct (8x)")
end

# ╔═╡ 8d6aa23c-521d-4348-a274-93c983f91465
md"""
Now this is looking even nicer! Both high and low frequencies have a similar amount of visual detail and the curve is looking smooth. 

But you may notice another annoyance: The fine noise around for example 2 kHz or 30 Hz. Certainly, this should have been smoothed away!

The reason for this remaining fine noise is that we used a rectangular window (`rect_window`) for our smoothing. This means all input points in a certain radius around an output point are weighted equally. If you now imagine moving the window just slightly, a point of different value could come in at the edge and change the value of the output, even though the window was just moved a bit. This is often also referred to as the rectangular window having large spectral leakage. Basically, it is not very good at removing fine (= high frequency) noise. 

The solution: Points at the edge must be weighed less than in the middle of the window. This is what we will do in the next section. 

## Using a Smoother Window

[`NonuniformResampling1D.jl`](https://github.com/Firionus/NonuniformResampling1D.jl) offers a selection of windows as an alternative to the rectangular window. A simple option is the Hann window, which is just a raised cosine function. Compared to the rectangular window, the Hann window is smooth and gives less weight to points that are further away:
"""

# ╔═╡ 19b0cc8b-c778-41cd-a140-e43da74f7e82
begin
	x = range(0, 1, length=100)
	plot(x, rect_window(.5).(x), label="rect_window(.5)", xlabel="position relative to next output point", title="Window Comparison")
	plot!(x, hann_window(1).(x), label="hann_window(1)")
end

# ╔═╡ e64cf445-613d-4917-b3d8-bfeff82d81f7
md"""
Note that for a similar smoothing effect, we have to use a wider width for the Hann window compared to the rectangular window. Here, we'll use twice the width. 

Applying this to our smoothing problem, the code looks like this:
"""

# ╔═╡ 55e4e82f-e0b9-47a3-bb85-102d35d36b09
oversampled_hann = nuresample(frange, magnitude, flog_oversampled, hann_window(1*oversampling))

# ╔═╡ afe5b231-6aad-4154-a346-26a9d7577a6a
begin
	plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Logarithmic Smoothing with Hann Window", label="raw", xlims=(20, 20e3), ylims=(40, 120), alpha=.3)
	plot!(flog_oversampled, oversampled_hann, label="log smooth 1/12 oct (8x, Hann 1.0)")
end

# ╔═╡ 7639a14c-449a-4611-806e-3405717debcf
md"""
As you can see, the Hann window strongly reduces the amount of high frequency noise in the smoothed output. But also note that it tends to make peaks look more soft and monotone, which may or may not be desirable. Lowering the width of the Hann window is an easy way to get back some more detail:
"""

# ╔═╡ d42f283f-afbf-4fab-bcb5-c97b5504a646
oversampled_hann2 = nuresample(frange, magnitude, flog_oversampled, hann_window(.7*oversampling))

# ╔═╡ b89209c2-f672-4099-9bbb-4f6022e85d81
begin
	plot(f, magnitude, xaxis=:log, xticks=[10, 100, 1e3, 10e3], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Sligthly Narrower Hann Window", label="raw", xlims=(20, 20e3), ylims=(40, 120), alpha=.3)
	plot!(flog_oversampled, oversampled_hann2, label="log smooth 1/12 oct (8x, Hann 0.7)")
end

# ╔═╡ 86ad9972-5280-42b2-85a6-434a12ee412c
md"""
> Comment: In the audio community, the amount of smoothing is usually communicated as just "1/12 octave smoothed", with no mention of the used window. This is because usually a rectangular window is assumed, so when you use a different window, you should document it. 
"""

# ╔═╡ 1edda514-9ada-4a9c-baef-2e8755a8bbe3
md"""
## Upsampling and Boundary Problems

For this part, let's focus on the super low frequency part of the measurement. 

> Comment: The measurement only contains noise below 30 Hz, but we will still use it for demonstration purposes. 

Let's plot the low frequencies around 1 Hz:
"""

# ╔═╡ dc5348fb-619e-4b9b-92cd-a10bf79db210
begin
	plot(f, magnitude, xaxis=:log, xticks=[0.1, 1, 10, 100], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Very Low Frequencies", label="raw", xlims=(.3, 20), ylims=(10, 90), alpha=1, legend=:topleft)
end

# ╔═╡ e2451a04-62f6-48bb-a7ee-a72e7f62c0fc
md"""
As you can see, on a log scale the density of points becomes quite low at low frequencies. Let's try to get a smoothed result:
"""

# ╔═╡ f4b54e68-5e9a-4975-a893-4c59ee74e1af
low_flog = octspace(.3, 20, resolution/oversampling)

# ╔═╡ 2fcae909-5596-4786-9336-749c14b573a1
low_smoothed = nuresample(frange, magnitude, low_flog, rect_window(.5*oversampling))

# ╔═╡ b9216b37-943d-4aca-96d4-15404cf8cf81
md"""
Oh no, an Error! It is telling us that there are not enough points at the beginning. Basically, we asked the resample algorithm to generate points around 0.3 Hz, where there is no input point in the window at all. 

This alone would not throw off the algorithm, since it will try to perform upsampling with a Lanczos3 algorithm and just interpolate points when there are not enough. However, interpolation with Lanczos3 takes 3 surrounding points to each side and these just don't exist in our case. In this case, `nuresample` will error instead of inventing values. 

The easiest way of dealing with this is to increase the lowest frequency for which we ask the resampling algorithm until the error goes away:
"""

# ╔═╡ e2ab1a01-b2a8-4686-bf09-9ba4d7285fb6
low_flog2 = octspace(.57, 20, resolution/oversampling)

# ╔═╡ cdebabe8-5bae-4995-b4f7-46b69e4ad640
low_smoothed2 = nuresample(frange, magnitude, low_flog2, rect_window(.5*oversampling))

# ╔═╡ 4aa18d3a-aaed-41df-aa8d-280e0e90be14
begin
	plot(f, magnitude, xaxis=:log, xticks=[0.1, 1, 10, 100], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Lanczos3 Interpolation", label="raw", xlims=(.3, 20), ylims=(10, 90), alpha=.3, legend=:topleft)
	plot!(low_flog2, low_smoothed2, label="log smooth 1/12 oct (8x)")
end

# ╔═╡ f24c21bc-5178-4c2c-b2d9-9d059ecdb915
md"""
And we get a beautiful result. You can see that the output has a similar level of detail everywhere, though it drops towards the very low frequencies, when the density of input points becomes less than 1/12 octaves. Still, the output is smooth even when there are little input points. This is because by default high quality Lanczos3 interpolation is performed before applying the smoothing. 

What other options are there to make a "not enough points" error go away? 

1. Avoid triggering the upsampling by changing the keyword argument `required_points_per_slice`. By default, this is set to 4, which means to each side of an output point there must be 4 input points in the window, otherwise 4 points will be interpolated.  
    But changing this value will not help much in our case, since the number of input points is so low and the output resolution is quite high, so upsampling is needed anyway to create at least one point in each window. 
2. Use a lower quality upsampling algorithm that needs less supporting points. For example, instead of Lanczos3 we could use linear interpolation:
"""

# ╔═╡ 285a45dd-2fd0-48aa-897a-57440e8c7001
low_smoothed_lin = nuresample(frange, magnitude, low_flog, rect_window(.5*oversampling), upsampling_function=tri_window()) # triangular window results in linear interpolation

# ╔═╡ e903849a-d56f-473b-82d6-5557baa416b9
begin
	plot(f, magnitude, xaxis=:log, xticks=[0.1, 1, 10, 100], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Linear Interpolation", label="raw", xlims=(.3, 20), ylims=(10, 90), alpha=.3, legend=:topleft)
	plot!(low_flog, low_smoothed_lin, label="log smooth 1/12 oct (8x)")
end

# ╔═╡ 35b07433-b177-48fe-a344-404c1b6bcd2d
md"""
Note that linear interpolation is applied before the normal 1/12 octave smoothing, so you can still see some slight differences between the linear interpolation in the raw data plot and the result of the log smoothing. 

Instead of using linear interpolation, you may visually prefer the look of Lanczos1 interpolation, which requires just as little surrounding points, but results in smoother output:
"""

# ╔═╡ 68de0a23-bdc7-4410-bef2-fe60b896e52f
low_smoothed4 = nuresample(frange, magnitude, low_flog, rect_window(.5*oversampling), upsampling_function=lanczos_window(1))

# ╔═╡ 51db438b-1e74-48b4-ad52-e025e183d254
begin
	plot(f, magnitude, xaxis=:log, xticks=[0.1, 1, 10, 100], minorticks=true, minorgrid=true, xlabel="f in Hz (log)", ylabel="dB", title="Lanczos1 Interpolation", label="raw", xlims=(.3, 20), ylims=(10, 90), alpha=.3, legend=:topleft)
	plot!(low_flog, low_smoothed4, label="log smooth 1/12 oct (8x)")
end

# ╔═╡ 160b6429-ad71-4503-9b05-1de81e898a2d
md"""
## Wrapping Up

I hope this demo of how to perform logarithmic smoothing with my packages [`NonuniformResampling1D`](https://github.com/Firionus/NonuniformResampling1D.jl) and [`NonlinearSequences`](https://github.com/Firionus/NonlinearSequences.jl) was useful to you. 

If you need other kinds of nonlinear smoothing, I hope it became clear that you can just replace `flog` with a differently spaced vector to get other kinds of nonlinear smoothing. `nuresample` is quite flexible there.

If you have any feedback regarding this demo or the used packages, feel free to open an issue in the corresponding GitHub repository. Even if it's just saying thank you or that you found a typo, an issue is a great way to reach out and I'd love to hear from you ❤️
"""

# ╔═╡ Cell order:
# ╟─e2ce0c30-0ce1-11ed-084b-4523774d8541
# ╟─ddfb9f83-a845-42d9-801a-d866528c3603
# ╠═c86650f9-4952-4064-80b6-e6d4c6253126
# ╠═dfb41db2-0b57-4de9-a44e-cd1f70842588
# ╟─cc004309-0ed9-465d-9f43-037ab1949a4e
# ╠═d8a3e8f1-83b4-425c-ace8-c48452c97c73
# ╟─2d0aba3c-b2e2-4ea1-8912-9e0d75532080
# ╠═1963ee15-c177-4e64-aa22-60fb783a1ec3
# ╟─df40186e-d36a-4779-be83-f709e102c43d
# ╠═35d2be48-81b5-4605-882a-684deb8675d3
# ╟─a4dfa7f8-b727-4936-bb24-63521168535a
# ╠═a44c0978-a733-4a90-b871-0c844586dcec
# ╟─80fab818-246f-4c69-bd43-faa23d4b3a16
# ╠═b255464e-366e-4329-b6c0-4d619d2be4ff
# ╟─9f6810e4-eadb-410b-850f-6bb9cda78800
# ╠═11f8da69-4a42-4f5c-a77c-11d94181c0a8
# ╟─5b1b6194-53f6-42b7-ac77-c691a373df3b
# ╠═dbfc2f92-d52a-4b04-b554-8d4c491d7549
# ╠═6def680a-5cff-474f-a556-eb311f118dba
# ╟─7956290b-3a19-40c2-8194-ff444d2b0e6a
# ╟─8a332d89-b63a-43a2-af02-907cfbd11505
# ╠═db47c349-4dc6-431a-8081-286506a674e3
# ╠═ebf1d7f5-ca3b-47fb-8fab-e49104a4fae5
# ╟─6d4e2046-4d64-4267-81a6-29df79763542
# ╠═7b56eec9-2717-49df-b39b-017b4d7a330c
# ╠═f67b1de1-1247-4ecf-be84-0bc1816bad34
# ╟─8d6aa23c-521d-4348-a274-93c983f91465
# ╠═19b0cc8b-c778-41cd-a140-e43da74f7e82
# ╟─e64cf445-613d-4917-b3d8-bfeff82d81f7
# ╠═55e4e82f-e0b9-47a3-bb85-102d35d36b09
# ╠═afe5b231-6aad-4154-a346-26a9d7577a6a
# ╟─7639a14c-449a-4611-806e-3405717debcf
# ╠═d42f283f-afbf-4fab-bcb5-c97b5504a646
# ╠═b89209c2-f672-4099-9bbb-4f6022e85d81
# ╟─86ad9972-5280-42b2-85a6-434a12ee412c
# ╟─1edda514-9ada-4a9c-baef-2e8755a8bbe3
# ╠═dc5348fb-619e-4b9b-92cd-a10bf79db210
# ╟─e2451a04-62f6-48bb-a7ee-a72e7f62c0fc
# ╠═f4b54e68-5e9a-4975-a893-4c59ee74e1af
# ╠═2fcae909-5596-4786-9336-749c14b573a1
# ╟─b9216b37-943d-4aca-96d4-15404cf8cf81
# ╠═e2ab1a01-b2a8-4686-bf09-9ba4d7285fb6
# ╠═cdebabe8-5bae-4995-b4f7-46b69e4ad640
# ╠═4aa18d3a-aaed-41df-aa8d-280e0e90be14
# ╟─f24c21bc-5178-4c2c-b2d9-9d059ecdb915
# ╠═285a45dd-2fd0-48aa-897a-57440e8c7001
# ╠═e903849a-d56f-473b-82d6-5557baa416b9
# ╟─35b07433-b177-48fe-a344-404c1b6bcd2d
# ╠═68de0a23-bdc7-4410-bef2-fe60b896e52f
# ╠═51db438b-1e74-48b4-ad52-e025e183d254
# ╟─160b6429-ad71-4503-9b05-1de81e898a2d
