function sph_hankel(z)
    # returns spherical hankel function 1st kind 1st order and its derivative
    h1  = exp.(1im*z)./z.*(1.0./(1im*z) .- 1.0)
    dh1 = exp.(1im*z)./z.*(-1im .+ 2.0./z .+ 2.0*1im./(z.^2))
    return h1,dh1
end

function exp_sph_hankel(z,zn)
    # returns spherical hankel function 1st kind 1st order and its derivative
    # multiplied with exp(-j*zn) in order to scale for large imaginary arguments
    exp_h1  = exp.(1im*(z .- zn))./z.*(1.0./(1im*z) .-1.0)
    exp_dh1 = exp.(1im*(z .- zn))./z.*(-1im .+ 2.0./z .+ 2.0*1im./(z.^2))
    return exp_h1,exp_dh1
end

function sphere_first_order(k,c,rho,a,up,rtp;S=-1,kv=Inf)

    # [p,v_r,v_theta,v_rA,v_thetaA,v_rV,v_thetaV] = sphere_fird_order(k,c,rho,a,up,rtp,S,kv);
    #
    # Sphere vibrating on its first order mode (up*cos.(theta))
    # See for example Lecture Note "Radiation of Sound", by  Finn Jacobsen
    # and Peter M. Juhl. See also Pierce.
    # If kv is included, the solution has viscous losses. The reference is
    # section 6.9 in 'Elements of Acoustics' by Temkin.
    #
    # Input:
    #   -k:        Wavenumber, m-1. Can be a vector.
    #   -c:        Speed of sound, m/s
    #   -a:        Sphere radius, m
    #   -rho:      Density of air, kg/m3
    #   -up:       Maximum normal velocity on the sphere, m/s
    #   -rtp:      Points where the sound pressure is calcualted. They are
    #              external to the sphere. It is a matrix with two columns,
    #              r and theta coordinates of the points, m
    #   -S:        sign convention (1, Karra exp(-jwt) or -1, Bruneau exp(jwt)).
    #              If not supplied, -1 is used (exp(jwt) convention).
    #   -kv:       if this parameter exists, the solution includes viscous
    #              losses with a viscous wavenumber kv. Can be a vector,
    #              following k.
    #
    # Ourput:
    #   -p:        Calculated sound pressure, Pa. Matrix with one row per
    #              point and one column per wavenumber.
    #   -v_r:      Total velocity in radial direction. Same size as p.
    #   -v_theta:  Total velocity in theta direction. Same size as p.
    #   -v_rA:     Acoustic (irrotational) velocity in radial direction. Same size as p.
    #   -v_thetaA: Acoustic (irrotational) velocity in theta direction. Same size as p.
    #   -v_rV:     Viscous (rotational) velocity in radial direction. Same size as p.
    #   -v_thetaV: Viscous (rotational) velocity in theta direction. Same size as p.
    #
    #  The last four outputs are empty if kv is not supplied (lossless calculation)

    # Peter Møller Juhl & Vicente Cutanda Henriquez 03-2012

    if any(rtp[:,1] .< a)
       error("There are points defined inside the sphere")
    end

    v_rA=[]
    v_thetaA=[]
    v_rV=[]
    v_thetaV=[]

    if kv != Inf
        # calculation including losses (uses exp(-jwt) convention)
        # tmp=zeros(1,length(kv))
        # tmp[1,:]=kv[1:end
        # kv=tmp; # kv becomes a row vector, in any case
        # impose exp(-jwt) convention to the viscous wavenumber
        kv=real(kv)+1im*abs.(imag(kv))

        Ka=kv*a # kv in VTconst correspond to K in Temkin
        ka=k*a # use perfect fluid wavenumber due to remarks above (6.9.1)
        # Note that here ka is wavenumber times radius - not acoustic wavenumber
        # as obtained from VTconst

        # Setup and solve (6.9.14-15)
        h1_ka, dh1_ka = sph_hankel(ka)
        h1_Ka, dh1_Ka = exp_sph_hankel(Ka,Ka) # scale to avoid very small coef's
        a11 = 3.0*1im*k.*dh1_ka
        a12 =-6.0*1im/a*h1_Ka
        a21 = 3.0*1im/a*h1_ka
        a22 =-3.0*1im/a*(Ka.*dh1_Ka+h1_Ka)

        A1=zeros(ComplexF64, size(k))
        exp_beta_B1=zeros(ComplexF64, size(k))
        exp_beta_B1=zeros(ComplexF64, size(k))
        AB=[a11 a12;
            a21 a22]\[up; up]
        A1=AB[1]
        exp_beta_B1=AB[2]
        A1=ones(ComplexF64,size(rtp,1))*A1
        exp_beta_B1=ones(ComplexF64,size(rtp,1))*exp_beta_B1

        # Find pressure potential (6.9.2)
        kr=rtp[:,1]*k
        h1_kr,dh1_kr=sph_hankel(kr)
        phi=1im*3.0*A1.*h1_kr.*(cos.(rtp[:,2])*ones(1,length(k))) # should become a matrix
        p=1im*rho*c*(ones(size(rtp,1),1)*k).*phi # (6.7.17)

        # Find radial velocity (6.9.11) with n=1
        Kr=rtp[:,1]*kv
        exp_h1_Kr,exp_dh1_Kr=exp_sph_hankel(Kr,ones(size(rtp,1),1)*Ka) #remember the scale on B1
        v_r=1im*3.0*(ones(size(rtp,1),1)*k).*(A1.*dh1_kr-2.0*exp_beta_B1./kr.*exp_h1_Kr).*(cos.(rtp[:,2])*ones(1,length(k)))

        v_rA=1im*3.0*(ones(size(rtp,1),1)*k).*(A1.*dh1_kr).*(cos.(rtp[:,2])*ones(1,length(k)))
        v_rV=1im*3.0*(ones(size(rtp,1),1)*k).*(-2.0*exp_beta_B1./kr.*exp_h1_Kr).*(cos.(rtp[:,2])*ones(1,length(k)))

        # Find theta component of velocity (6.9.13) with n=1
        v_theta=1im*3.0./(rtp[:,1]*ones(1,length(k))).*(A1.*h1_kr-exp_beta_B1.*(Kr.*exp_dh1_Kr+exp_h1_Kr)).*(-sin.(rtp[:,2])*ones(1,length(k)))

        v_thetaA=1im*3.0./(rtp[:,1]*ones(1,length(k))).*(A1.*h1_kr).*(sin.(rtp[:,2])*ones(1,length(k)))
        v_thetaV=1im*3.0./(rtp[:,1]*ones(1,length(k))).*(-exp_beta_B1.*(Kr.*exp_dh1_Kr+exp_h1_Kr)).*(sin.(rtp[:,2])*ones(1,length(k)))


        if S==-1 # Change to exp(jwt) convention if required
            p=conj(p)
            v_r=conj(v_r); v_theta=conj(v_theta)
            v_rA=conj(v_rA); v_thetaA=conj(v_thetaA)
            v_rV=conj(v_rV); v_thetaV=conj(v_thetaV)
        end

    else
        # calculation without losses  (uses exp(jwt) convention)
        A1=-ones(size(rtp,1),1)*(rho*c*a*up*k.*exp(1i*k*a)./(1 - 2.0./(k*a).^2 + 2.0./(j*k*a)))  # (5.38)

        p=(A1.*(1 + 1.0./(1i*rtp[:,1]*k) )).*(exp.(-1i*rtp[:,1]*k)./(rtp[:,1]*k)).*
            (cos.(rtp[:,2])*ones(1,length(k))) # FJ (jwt convention)
        v_r=-A1/(rho*c).*(1-2.0./(rtp[:,1]*k).^2+2.0./(1im*rtp[:,1]*k)).*(exp.(-1im.*
            (rtp[:,1]*k))./(rtp[:,1]*k)).*(cos.(rtp[:,2])*ones(1,length(k))) #(5.35)
        v_theta=(1im.*A1)/(rho*c).*(1+1.0./(1im*rtp[:,1]*k)).*
            (exp.(-1im*rtp[:,1]*k)./((rtp[:,1].^2)*k)).*(sin.(rtp[:,2])*ones(1,length(k))) # (5.36)

        if S==1 # Change to exp(-jwt) convention if required
            p=conj(p)
            v_r=conj(v_r); v_theta=conj(v_theta)
        end

    end

    return p, v_r, v_theta, v_rA, v_thetaA, v_rV, v_thetaV

end
