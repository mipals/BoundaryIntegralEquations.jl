pressure(mn,ka) = (mn == 0 ? atan(-besselj(1,ka)/bessely(1,ka)) :
            atan((besselj(mn-1,ka)-besselj(mn+1,ka))/(bessely(mn+1,ka)-bessely(mn-1,ka))))
function plane_wave_scattering_circle(ϕ,ka,nterm=10)
    int_incident,int_scattering   = (ones(length(ϕ))*besselj(0,ka), zeros(length(ϕ)))
    for m=1:nterm-1; int_incident = int_incident + 2*im .^(m) * besselj(m,ka) .* cos.(m*(ϕ)) end
    for m=0:nterm-1;
        int_scattering = int_scattering - (m == 0 ? 1.0 : 2.0)*im .^(m+1.0) * exp(-im*pressure(m,ka)) *
                    sin.(pressure(m,ka)) * (besselj(m,ka)+im*bessely(m,ka)) .* cos.(m*(ϕ))
    end
    return int_incident + int_scattering
end
