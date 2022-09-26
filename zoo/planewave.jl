struct PlaneWave{N,T}
    k::SVector{N,T}
    ω::T
end

wavevector(x::PlaneWave) = x.k
frequency(x::PlaneWave) = x.ω
