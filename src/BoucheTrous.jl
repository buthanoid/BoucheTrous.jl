module BoucheTrous

using SimpleAlgebra, AbstractFFTs,OptimPackNextGen, Zygote, EasyFITS

export BoucheTrous!, BoucheTrous_IRDIS

function BoucheTrous!(
    img::AbstractArray{T,N}, badpix::BitArray{N}, pupil_radius::AbstractFloat
) where {T,N}

    sz = size(img)
    pupil = sqrt.(rfftfreq(sz[1]).^2 .+ (fftfreq(sz[2]).^2)') .< pupil_radius

    F = LinOpDFT(T,sz[1:2])
    #pupil[(1:n for n ∈ outputsize(F))...]
    D = LinOpDiag(T,.!pupil)
    C = CostL2(T,size(pupil),0)

    Threads.@threads for n=1:size(img,3)
        data = img[:,:,n]
        S = LinOpSelect(T,badpix[:,:,n])
        H = C * D *F * (S' + data)
        cost(x) = H*x
        p = vmlmb( cost, S*data; lower=0., autodiff=true, xtol=(0,0.05), verb=false)
        img[:,:,n] .+=  S'*p
    end
end

# diameter [m], λ [m]
function compute_pupil_radius(diameter, λ, pixscal) ::Float64
    RAD2MAS = 1 / π *180 *60*60*1000
    pupil_radius = pixscal / RAD2MAS  * diameter / λ
end

# [nm]
const BANDS_DICT_IRDIS = Dict{String,Vector{Float64}}(
    # DBI
    "DB_Y23" => [1.022, 1.076], "DB_J23" => [1.190, 1.273], "DB_H23" => [1.593, 1.667],
    "DB_H32" => [1.667, 1.733], "DB_H34" => [1.667, 1.733], "DB_K12" => [2.110, 2.251],
    # DBI + ND filter
    "DB_NDH23" => [1.593, 1.667], "DB_NDH32" => [1.667, 1.733], #TODO: is it correct ?
    # Broadband
    "BB_Y" => [1.043, 1.043], "BB_J"  => [1.245, 1.245],
    "BB_H" => [1.625, 1.625], "BB_Ks" => [2.182, 2.182],
    # DP Broadband #TODO: are all DP present and correct ?
    "DP_0_BB_Y" => [1.043, 1.043], "DP_0_BB_J"  => [1.245, 1.245],
    "DP_0_BB_H" => [1.625, 1.625], "DP_0_BB_Ks" => [2.182, 2.182],
    # Narrowband
    "NB_HeI"    => [1.085, 1.085], "NB_ContJ" => [1.213, 1.213], "NB_PaB"    => [1.283, 1.283],
    "NB_ContH"  => [1.573, 1.573], "NB_FeII"  => [1.642, 1.642], "NB_ContK1" => [2.091, 2.091],
    "NB_H2"     => [2.124, 2.124], "NB_BrG"   => [2.170, 2.170], "NB_CO"     => [2.290, 2.290],
    "NB_ContK2" => [2.266, 2.266],
    # Spectroscopy
    "S_LR" => [1.644,  1.644], "S_MR" => [1.3845, 1.3845])
    

const IGNORED_PIXELS_IRDIS = begin
    vpm = falses(2048, 1024)
    vpm[   1:50   ,     :    ] .= true
    vpm[    :     ,    1:18  ] .= true
    vpm[1962:2048 ,     :    ] .= true
    vpm[    :     , 1014:1024] .= true
    vpm[ 935:1080 ,     :    ] .= true
    vpm
end

"""
    BoucheTrous_IRDIS(::FitsFile; write=false, kwds...) -> Array

Compute `BoucheTrous!` for the left and right image of each frame of the given IRDIS file, and
optionaly write the result in a new FITS file.

# Keyword parameters
- `extdata    ::Union{Integer,AbstractString} = 1`: the extension of the HDU containing the data.
  Default is `1` which means the primary HDU.
- `extweights ::Union{Integer,AbstractString} = "weights"`: the extension of the HDU containing
  the weights. This gives the bad pixels to replace, which are the pixels with a weight of `0`.
- `write ::Bool = false`: when `true`, will write the result in a new file, with the data HDU
  replaced with the "bouched" data and with a keyword in the header.
- `dst_filepath ::String =""`: only used when `write` is `true`. When empty, the new file will be
  created in the same folder of the src file, with a prefix. When non empty, the given `String`
  will be used as filepath.
"""
function BoucheTrous_IRDIS(
    fitsfile::FitsFile
    ; extdata=1, extweights="weights", write::Bool=false, dst_filepath::String=""
) ::Array

    hdudata = fitsfile[extdata]
    hduweights = fitsfile[extweights]
    
    diameter = 8.2 # diametre des UT [m]
    
    band = hdudata["ESO INS COMB IFLT"].string
    (λl, λr) = BANDS_DICT_IRDIS[band] .* 1e-6 # [m]
    
    
    pixscal = hdudata["PIXSCAL"].float

    pupil_radius_l = compute_pupil_radius(diameter, λl, pixscal)
    pupil_radius_r = compute_pupil_radius(diameter, λr, pixscal)

    data = read(hdudata)
    badpixels ::BitArray = (read(hduweights) .== 0)
    
    # we do not care about some of the bad pixels
    badpixels[IGNORED_PIXELS_IRDIS,:] .= false
    
    data_l      = view(data, 1:1024,    1:1024, :)
    data_r      = view(data, 1025:2048, 1:1024, :)
    badpixels_l = badpixels[1:1024,    1:1024, :]
    badpixels_r = badpixels[1025:2048, 1:1024, :]
    
    BoucheTrous!(data_l, badpixels_l, pupil_radius_l)
    BoucheTrous!(data_r, badpixels_r, pupil_radius_r)

    if write
        if isempty(dst_filepath)
            src_filepath = fitsfile.path
            dst_filepath = joinpath(dirname(src_filepath), "bouched_" * basename(src_filepath))
        end
        hdr = FitsHeader(hdudata)
        hdr["BOUCHTRU"] = (true, "bad pixels replaced")
        modify_one_hdu(fitsfile, dst_filepath, extdata, hdr, data)
    end

    return data
end

"""
    BoucheTrous_IRDIS(fitspath::AbstractString; kwds...)

Alternative version with a FITS file path instead of a `FitsFile` handler.
"""
function BoucheTrous_IRDIS(fitspath::AbstractString; kwds...)
    FitsFile(path -> BoucheTrous_IRDIS(path; kwds...), fitspath)
end

"""
    modify_one_hdu(src_fitsfile::FitsFile, dst_filepath::AbstractString, \
        ext::Union{Integer,AbstractString}, new_hdr::FitsHeader, new_data::Array) -> Bool

Write a copy of `src_fitsfile` at `dst_filepath` but with hdu `ext` replaced with `new_hdr` and \
`new_data`.

The destination cannot be the source. The HDU to modify must exist and be a `FitsImageHDU`.
"""
function modify_one_hdu(
    src_fitsfile::FitsFile, dst_filepath::AbstractString, ext::Union{Integer,AbstractString},
    new_hdr::FitsHeader, new_data::Array
) ::Nothing

    src_fitsfile.path == dst_filepath && throw(ArgumentError("dst path cannot equal src path"))
    
    hdu_to_modify ::FitsImageHDU = src_fitsfile[ext]
    
    modified = false
    
    FitsFile(dst_filepath, "w!") do dst_fitsfile
        for hdu in src_fitsfile
            if !modified && hdu == hdu_to_modify
                write(dst_fitsfile, new_hdr, new_data)
                modified = true
            else
                write(dst_fitsfile, FitsHeader(hdu), read(hdu))
            end
        end
    end
    
    nothing
end

end # module BoucheTrous


